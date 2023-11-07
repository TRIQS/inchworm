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
  read_json_parameters("/Users/yangyu/src/inchworm/apps/parameters.json", debug, cp, n_site, epsilon, theta, n_bath, n_spin, U, mu, t, tau_max,
                       tau_split, n_GK, bond_dim, sweep_bound, order_list, tci_prrlu, error_bound, bl_index, subspace_index);
  auto [vi, wi] = select_quadrature_GK(n_GK, 0, 1);

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

  std::ofstream outfile("output.txt");
  // TCI
  int n_phi                                 = std::accumulate(block_shape.begin(), block_shape.end(), 0);
  std::vector<double> integral_order_list   = {};
  std::vector<double> calculation_time_list = {};
  for (int order : order_list) {
    auto start_time = std::chrono::high_resolution_clock::now();
    int n           = 2 * order; // number of tau's
    std::vector<int> pivot1(n, 0);
    std::vector<int> range(n);
    std::iota(range.begin(), range.end(), 0);
    auto phi_pair_list      = get_all_phi_crossing(range);      //gives all possible phi
    auto iota_pair_list     = get_all_iota(block_shape, order); //gives all possible iota
    auto iota1              = pivot1;
    double integral_sum_phi = 0.0;
    // for (auto [phi_d_list, phi_d_dag_list] : phi_pair_list) { // phi start
    auto [phi_d_list, phi_d_dag_list] = phi_pair_list[0];
    double integral_sum_n_left        = 0.0;
    // for (int n_left = 1; n_left < n; n_left++) { // n_left start
    double integral_sum_iota = 0.0;
    int n_left               = order;
    // for (auto [iota_d_list, iota_d_dag_list] : iota_pair_list) { // iota start
    long count = 0;
    std::vector<double> vs1;
    for (int i = 0; i < pivot1.size(); i++) { vs1.push_back(vi[pivot1[i]]); }
    // auto get_u_tau_max_00 = [&u_tau_max_zeroth_order, &tau_split, &tau_max, &all_d_ops, &all_d_dag_ops, &block_shape, &cp, &Delta_tau = Delta_tau,
    //                          &ad_imp = ad_imp, &u_interpolator, &count, &phi_d_list = phi_d_list, &phi_d_dag_list = phi_d_dag_list,
    //                          &iota_d_list = iota_d_list, &iota_d_dag_list = iota_d_dag_list, &n_left, &bl_index, &subspace_index,
    //                          &debug](const std::vector<double> &vs) {
    //   std::vector<double> vs_left(vs.begin(), vs.begin() + n_left);
    //   std::vector<double> vs_right(vs.begin() + n_left, vs.end());
    //   std::vector<double> taus_left  = changeVariable(vs_left, tau_split, 0.0);
    //   std::vector<double> taus_right = changeVariable(vs_right, tau_max, tau_split);
    //   auto taus(taus_left);
    //   taus.insert(taus.end(), taus_right.begin(), taus_right.end());
    //   double integrand = evaluate_u_tau_max(u_tau_max_zeroth_order, tau_split, tau_max, all_d_ops, all_d_dag_ops, block_shape, cp, Delta_tau,
    //                                         ad_imp, u_interpolator, getElements(phi_d_list, taus), getElements(phi_d_dag_list, taus), iota_d_list,
    //                                         iota_d_dag_list, bl_index, subspace_index);
    //   count++;
    //   double j = jacobian(taus_left, tau_split, 0.0) * jacobian(taus_right, tau_max, tau_split);
    //   return integrand * j;
    // };
    auto vs               = vs1;
    auto get_u_tau_max_00 = [&u_tau_max_zeroth_order, &tau_split, &tau_max, &all_d_ops, &all_d_dag_ops, &block_shape, &cp, &Delta_tau = Delta_tau,
                             &ad_imp = ad_imp, &u_interpolator, &count, &phi_d_list = phi_d_list, &phi_d_dag_list = phi_d_dag_list, &n_left,
                             &bl_index, &subspace_index, &debug, &vs, &n](const std::vector<int> &iota_list) {
      std::vector<double> vs_left(vs.begin(), vs.begin() + n_left);
      std::vector<double> vs_right(vs.begin() + n_left, vs.end());
      std::vector<double> taus_left  = changeVariable(vs_left, tau_split, 0.0);
      std::vector<double> taus_right = changeVariable(vs_right, tau_max, tau_split);
      auto taus(taus_left);
      std::vector<int> iota_d_list(iota_list.begin(), iota_list.begin() + n / 2);
      std::vector<int> iota_d_dag_list(iota_list.begin() + n / 2, iota_list.end());
      // std::cout << "iota_d_list: ";
      // print_vector(iota_d_list);
      // std::cout << "iota_d_dag_list: ";
      // print_vector(iota_d_dag_list);
      taus.insert(taus.end(), taus_right.begin(), taus_right.end());
      double integrand =
         evaluate_u_tau_max(u_tau_max_zeroth_order, tau_split, tau_max, all_d_ops, all_d_dag_ops, block_shape, cp, Delta_tau, ad_imp, u_interpolator,
                            getElements(phi_d_list, taus), getElements(phi_d_dag_list, taus), iota_d_list, iota_d_dag_list, bl_index, subspace_index);
      count++;
      double j = jacobian(taus_left, tau_split, 0.0) * jacobian(taus_right, tau_max, tau_split);
      return integrand * j;
    };
    double u_tau_max_00_vs1 = 0;
    if (debug) {
      std::vector<double> vs1_left(vs1.begin(), vs1.begin() + n_left);
      std::vector<double> vs1_right(vs1.begin() + n_left, vs1.end());
      std::vector<double> taus1_left  = changeVariable(vs1_left, tau_split, 0.0);
      std::vector<double> taus1_right = changeVariable(vs1_right, tau_max, tau_split);
      auto taus(taus1_left);
      taus.insert(taus.end(), taus1_right.begin(), taus1_right.end());
      // std::cout << "iota_d_list: ";
      // print_vector(iota_d_list);
      // std::cout << "iota_d_dag_list: ";
      // print_vector(iota_d_dag_list);
      std::cout << "tau_d_list: ";
      print_vector(getElements(phi_d_list, taus));
      std::cout << "tau_d_dag_list: ";
      print_vector(getElements(phi_d_dag_list, taus));
      u_tau_max_00_vs1 = get_u_tau_max_00(iota1);
      // u_tau_max_00_vs1 = get_u_tau_max_00(vs1);
      std::cout << "get_u_tau_max(pivot1): " << u_tau_max_00_vs1 << "\n" << std::endl;
      // for (auto iota_d : iota_d_list) { outfile << iota_d; }
      // for (auto iota_d_dag : iota_d_dag_list) { outfile << iota_d_dag; }
      std::cout << "iota1: ";
      print_vector(iota1);
      outfile << " ";
      outfile << std::fixed << std::setprecision(20) << u_tau_max_00_vs1 << "\n";
    } else {
      u_tau_max_00_vs1 = get_u_tau_max_00(iota1);
    }
    if (u_tau_max_00_vs1 == 0) { continue; }

    double current_integral{0};
    double previous_integral{0};
    if (debug) { std::cout << "iteration nEval LastSweepPivotError integral\n"; }
    wi = std::vector(n_phi, 1.0);
    std::vector<int> iotai(n_phi);
    std::iota(iotai.begin(), iotai.end(), 0);
    if (tci_prrlu) {
      auto ci = xfac::CTensorCI2<double, int>(get_u_tau_max_00, std::vector(n, iotai), {.bond_dim = bond_dim, .pivot1 = pivot1});
      for (int i = 0; i < sweep_bound; i++) {
        ci.iterate();
        ci.makeCanonical();
        current_integral = ci.tt.sum(std::vector(n, wi));
        if (debug) { std::cout << i << " " << count << " " << ci.pivotError[ci.pivotError.size() - 1] << " " << current_integral << std::endl; }
        if (std::abs(current_integral - previous_integral) < error_bound && i > 1) { break; }
        previous_integral = current_integral;
      }
      std::cout <<"True error:"<< ci.tt.trueError(get_u_tau_max_00) << std::endl;
    } else {
      auto ci = xfac::CTensorCI<double, int>(get_u_tau_max_00, std::vector(n, iotai), {.pivot1 = pivot1});
      for (int i = 0; i < sweep_bound; i++) {
        ci.iterate();
        current_integral = ci.sumWeighted(std::vector(n, wi));
        if (debug) { std::cout << i << " " << count << " " << ci.pivotError[ci.pivotError.size() - 1] << " " << current_integral << std::endl; }
        if (std::abs(current_integral - previous_integral) < error_bound && i > 1) { break; }
        previous_integral = current_integral;
      }
    }
    if (debug) { std::cout << std::endl; }
    integral_sum_iota += current_integral;
    // } // iota end
    integral_sum_n_left += integral_sum_iota;
    // } // n_left end
    integral_sum_phi += integral_sum_n_left;
    // } // phi end
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
