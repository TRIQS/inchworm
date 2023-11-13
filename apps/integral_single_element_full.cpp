#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>

#include <xfac/grid.h>
#include <xfac/tensor/tensor_ci.h>
#include <xfac/tensor/tensor_ci_2.h>
#include <xfac/tensor/tensor_train.h>
#include <inchworm/diagram/diagram.hpp>
#include <inchworm/diagram/inclusion_exclusion.hpp>
#include <inchworm/diagram/print.hpp>
#include <inchworm/atom_diag.hpp>
#include <inchworm/u_frame.hpp>
#include <inchworm/impurity_product.hpp>
#include <inchworm/util.hpp>
#include <inchworm/interpolator.hpp>
#include "./hubbard.hpp"
#include "./integral_util.hpp"
using namespace inchworm;

int main() {

  // parameters
  bool debug = false;
  constr_params_t cp{};
  mat_t theta;
  vec_t epsilon;
  int n_site{};
  int n_bath{};
  int n_spin{};
  int n_GK{};
  int bond_dim{};
  int sweep_bound{};
  double U{};
  double mu{};
  double t{};
  double tau_max{};
  double tau_split{};
  std::vector<int> order_list = {};
  bool tci_prrlu{};
  double error_bound{};
  int bl_index{};
  int subspace_index{};
  read_json_parameters("../../apps/parameters.json", debug, cp, n_site, epsilon, theta, n_bath, n_spin, U, mu, t, tau_max,
                       tau_split, n_GK, bond_dim, sweep_bound, order_list, tci_prrlu, error_bound, bl_index, subspace_index);

  // prepare input
  auto [Delta_tau, ad_imp, u_tau, G_tau]        = test_setup(n_site, n_bath, n_spin, U, mu, t, cp, theta, epsilon);
  auto u_interpolator                           = interpolator_t<scalar_t>(u_tau, u_tau[0].mesh().size());
  frame_t u_tau_max_zeroth_order                = u_interpolator(tau_max - tau_split) * u_interpolator(tau_split); //oder 0 result
  long n_bl                                     = cp.gf_struct.size();
  std::vector<int> block_shape                  = {};
  std::vector<std::vector<fop_t>> all_d_ops     = {};
  std::vector<std::vector<fop_t>> all_d_dag_ops = {};
  all_d_ops.resize(n_bl, {});
  all_d_dag_ops.resize(n_bl, {});
  for (auto [bl, bl_pair] : enumerate(cp.gf_struct)) {
    auto [bl_name, bl_size] = bl_pair;
    block_shape.push_back(bl_size);
  }
  auto fops = fundamental_operator_set{cp.gf_struct};
  for (auto [bl, bl_pair] : enumerate(cp.gf_struct)) {
    auto [bl_name, bl_size] = bl_pair;
    all_d_ops[bl].clear();
    all_d_dag_ops[bl].clear();
    for (auto idx : range(bl_size)) {
      all_d_ops[bl].emplace_back(0.0, false, fops[{bl_name, idx}], bl, idx);
      all_d_dag_ops[bl].emplace_back(0.0, true, fops[{bl_name, idx}], bl, idx);
    }
  }
  std::cout << "Delta_tau shape:" << std::endl;
  print_block_shape(Delta_tau);
  std::cout << "G_tau shape:" << std::endl;
  print_block_shape(G_tau);
  std::cout << "u_tau shape:" << std::endl;
  print_block_shape(u_tau);
  std::cout << std::setprecision(17) << std::endl;

  // TCI
  std::vector<double> integral_order_list   = {};
  std::vector<double> calculation_time_list = {};
  auto [vi, wi_v]                           = select_quadrature_GK(n_GK, 0, 1);
  int n_phi                                 = std::accumulate(block_shape.begin(), block_shape.end(), 0);
  std::vector<double> iotai(n_phi);
  std::iota(iotai.begin(), iotai.end(), 0);
  auto wi_iota = std::vector(n_phi, 1.0);
  for (int order : order_list) {
    auto start_time = std::chrono::high_resolution_clock::now();
    int n           = 2 * order;     // number of tau's
    std::vector<int> v_pivot1(n, 0); // for tau only; pivots for iota are set later
    std::vector<int> range(n);
    std::iota(range.begin(), range.end(), 0);
    auto phi_pair_list = get_all_phi(range); //gives all possible phi

    std::vector<int> iota_pivots(n, 0);
    std::vector<int> iota_pivots_range(n_phi);
    std::iota(iota_pivots_range.begin(), iota_pivots_range.end(), 0);
    std::vector<std::vector<int>> all_iota_pivots{};
    generate_combinations(iota_pivots_range, iota_pivots, 0, all_iota_pivots);
    // auto iota_pair_list     = get_all_iota(block_shape, order); //gives all possible iota
    double integral_sum_phi = 0.0;
    for (auto [phi_d_list, phi_d_dag_list] : phi_pair_list) {
      double integral_sum_n_left = 0.0;
      for (int n_left = 1; n_left < n; n_left++) {
        double integral_sum_iota = 0.0;
        // for (auto [iota_d_list, iota_d_dag_list] : iota_pair_list) {
        long count                 = 0;
        auto get_u_tau_max_element = [&u_tau_max_zeroth_order, &tau_split, &tau_max, &all_d_ops, &all_d_dag_ops, &block_shape, &cp,
                                      &Delta_tau = Delta_tau, &ad_imp = ad_imp, &u_interpolator, &count, &phi_d_list = phi_d_list,
                                      &phi_d_dag_list = phi_d_dag_list, &n_left, &bl_index, &subspace_index,
                                      &debug](const std::vector<double> &v_iota_s) {
          std::vector<double> vs{};
          std::vector<double> iotas{};
          vs.reserve(v_iota_s.size() / 2);
          iotas.reserve(v_iota_s.size() / 2);
          for (int i = 0; i < v_iota_s.size(); i++) {
            if (i % 2 == 0) {
              iotas.push_back(v_iota_s[i]);
            } else {
              vs.push_back(v_iota_s[i]);
            }
          }
          // std::cout << "vs: ";
          // print_vector(vs);
          // std::cout << "iotas: ";
          // print_vector(iotas);
          // for(int i = 0 ; i < iotas.size(); i++){
          //   if(iotas[i]!=0){
          //     std::cout << "iota is not zero" << std::endl;
          //   }
          // }
          // int mid_iota = iotas.size() / 2;
          // std::vector<double> iota_d_list(iotas.begin(), iotas.begin() + mid_iota);
          // std::vector<double> iota_d_dag_list(iotas.begin() + mid_iota, iotas.end());
          std::vector<double> iota_d_list     = getElements(phi_d_list, iotas);
          std::vector<double> iota_d_dag_list = getElements(phi_d_dag_list, iotas);
          std::vector<double> vs_left(vs.begin(), vs.begin() + n_left);
          std::vector<double> vs_right(vs.begin() + n_left, vs.end());
          std::vector<double> taus_left  = changeVariable(vs_left, tau_split, 0.0);
          std::vector<double> taus_right = changeVariable(vs_right, tau_max, tau_split);
          auto taus(taus_left);
          taus.insert(taus.end(), taus_right.begin(), taus_right.end());
          double integrand = evaluate_u_tau_max(u_tau_max_zeroth_order, tau_split, tau_max, all_d_ops, all_d_dag_ops, block_shape, cp, Delta_tau,
                                                ad_imp, u_interpolator, getElements(phi_d_list, taus), getElements(phi_d_dag_list, taus), iota_d_list,
                                                iota_d_dag_list, bl_index, subspace_index);
          count++;
          double j = jacobian(taus_left, tau_split, 0.0) * jacobian(taus_right, tau_max, tau_split);
          return integrand * j;
        };

        std::vector<double> v_iota_s1;
        for (int i = 0; i < v_pivot1.size(); i++) { v_iota_s1.push_back(vi[v_pivot1[i]]); }
        double u_tau_max_element_vs1 = 0;

        // set pivot for iota
        int iota_pivot_index = 0;
        for (auto iota_pivot1 : all_iota_pivots) {
          std::vector<double> v_iota_s1_temp{};
          v_iota_s1_temp.reserve(v_iota_s1.size() + iota_pivot1.size());
          for (int i = 0; i < iota_pivot1.size(); i++) {
            v_iota_s1_temp.push_back(iota_pivot1[i]);
            v_iota_s1_temp.push_back(v_iota_s1[i]);
          }
          u_tau_max_element_vs1 = get_u_tau_max_element(v_iota_s1_temp);
          if (u_tau_max_element_vs1 != 0) {
            v_iota_s1 = v_iota_s1_temp;
            break;
          }
          iota_pivot_index++;
        }
        if (iota_pivot_index == all_iota_pivots.size()) { continue; }

        std::vector<int> pivot1{};
        auto pivot1_to_append = all_iota_pivots[iota_pivot_index];
        pivot1.reserve(v_pivot1.size() + pivot1_to_append.size());
        for (int i = 0; i < pivot1_to_append.size(); i++) {
          pivot1.push_back(pivot1_to_append[i]);
          pivot1.push_back(v_pivot1[i]);
        }
        std::cout << "pivot1: ";
        print_vector(pivot1);
        std::cout << "v_iota_s1: ";
        print_vector(v_iota_s1);

        if (debug) {
          std::vector<double> vs1{};
          std::vector<double> iotas1{};
          vs1.reserve(v_iota_s1.size() / 2);
          iotas1.reserve(v_iota_s1.size() / 2);
          for (int i = 0; i < v_iota_s1.size(); i++) {
            if (i % 2 == 0) {
              iotas1.push_back(v_iota_s1[i]);
            } else {
              vs1.push_back(v_iota_s1[i]);
            }
          }
          int mid_iota1 = iotas1.size() / 2;
          std::vector<double> iota_d_list1(iotas1.begin(), iotas1.begin() + mid_iota1);
          std::vector<double> iota_d_dag_list1(iotas1.begin() + mid_iota1, iotas1.end());
          std::vector<double> vs1_left(vs1.begin(), vs1.begin() + n_left);
          std::vector<double> vs1_right(vs1.begin() + n_left, vs1.end());
          std::vector<double> taus1_left  = changeVariable(vs1_left, tau_split, 0.0);
          std::vector<double> taus1_right = changeVariable(vs1_right, tau_max, tau_split);
          auto taus(taus1_left);
          taus.insert(taus.end(), taus1_right.begin(), taus1_right.end());
          std::cout << "iota_d_list: ";
          print_vector(iota_d_list1);
          std::cout << "iota_d_dag_list: ";
          print_vector(iota_d_dag_list1);
          std::cout << "tau_d_list: ";
          print_vector(getElements(phi_d_list, taus));
          std::cout << "tau_d_dag_list: ";
          print_vector(getElements(phi_d_dag_list, taus));
          std::cout << "get_u_tau_max_element(pivot1): " << u_tau_max_element_vs1 << "\n" << std::endl;
        }
        if (u_tau_max_element_vs1 == 0) { continue; }

        double current_integral{0};
        double previous_integral{0};
        double integral_element{0};
        if (debug) { std::cout << "iteration nEval LastSweepPivotError integral\n"; }
        std::vector<std::vector<double>> input{};
        input.reserve(2 * n);
        for (int i = 0; i < n; i++) {
          input.push_back(iotai);
          input.push_back(vi);
        }
        std::vector<std::vector<double>> weight{};
        weight.reserve(2 * n);
        for (int i = 0; i < n; i++) {
          weight.push_back(wi_iota);
          weight.push_back(wi_v);
        }
        // std::cout << "input: " << std::endl;
        // for (auto v : input) { print_vector(v); }
        // std::cout << "weight: " << std::endl;
        // for (auto v : weight) { print_vector(v); }
        if (tci_prrlu) {
          auto ci = xfac::CTensorCI2<double, double>(get_u_tau_max_element, input, {.bond_dim = bond_dim, .pivot1 = pivot1});
          for (int i = 0; i < sweep_bound; i++) {
            ci.iterate();
            // ci.makeCanonical();
            current_integral = ci.tt.sum(weight);
            if (debug) { std::cout << i << " " << count << " " << ci.pivotError[ci.pivotError.size() - 1] << " " << current_integral << std::endl; }
            if (std::abs(current_integral - previous_integral) < error_bound && i > 1) { break; }
            previous_integral = current_integral;
          }
          integral_element = current_integral;
          if (debug) { print_rank(ci.tt); }
        } else {
          auto ci = xfac::CTensorCI<double, double>(get_u_tau_max_element, input, {.pivot1 = pivot1});
          for (int i = 0; i < sweep_bound; i++) {
            ci.iterate();
            current_integral = ci.sumWeighted(weight);
            if (debug) { std::cout << i << " " << count << " " << ci.pivotError[ci.pivotError.size() - 1] << " " << current_integral << std::endl; }
            if (std::abs(current_integral - previous_integral) < error_bound && i > 5) { break; }
            previous_integral = current_integral;
          }
          integral_element = current_integral;
          if (debug) {
            std::cout << "rank:" << std::endl;
            print_vector(ci.rank());
          }
        }
        if (debug) { std::cout << std::endl; }
        integral_sum_iota += integral_element;
        // }
        integral_sum_n_left += integral_sum_iota;
      }
      integral_sum_phi += integral_sum_n_left;
    }
    auto end_time            = std::chrono::high_resolution_clock::now();
    auto duration            = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    auto duration_in_seconds = static_cast<double>(duration) / 1e6;
    calculation_time_list.push_back(duration_in_seconds);
    integral_order_list.push_back(integral_sum_phi);
  }

  // print results
  int i = subspace_index / u_tau[bl_index].target_shape()[0];
  int j = subspace_index % u_tau[bl_index].target_shape()[0];
  std::cout << "u_tau_max exact: " << std::setw(10) << u_interpolator(tau_max)[bl_index](i, j) << std::endl;
  std::cout << std::left << std::setw(10) << "order" << std::setw(30) << "value" << std::setw(30) << "time(s)" << std::endl;
  std::cout << std::setw(10) << "0" << std::setw(30) << u_tau_max_zeroth_order[bl_index](i, j) << std::endl;
  for (int i = 0; i < order_list.size(); i++) {
    std::cout << std::setw(10) << order_list[i] << std::setw(30) << integral_order_list[i] << std::setw(30) << calculation_time_list[i] << std::endl;
  }
  double sum_value = u_tau_max_zeroth_order[bl_index](i, j) + std::accumulate(integral_order_list.begin(), integral_order_list.end(), 0.0);
  double sum_time  = std::accumulate(calculation_time_list.begin(), calculation_time_list.end(), 0.0);
  std::cout << std::setw(10) << "sum:" << std::setw(30) << sum_value << std::setw(10) << sum_time << std::endl;

  return 0;
}
