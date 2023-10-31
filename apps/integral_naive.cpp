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
#include "./integral_util_naive.hpp"
using namespace inchworm;

int main() {

  //parameters for the model
  constr_params_t cp;
  cp.beta          = 2.0;
  cp.gf_struct     = {{"up", 1}};
  cp.n_tau_green   = 5;
  cp.n_tau_inch    = 10001;
  cp.n_tau         = 10001;
  mat_t theta      = {{1.0}};
  vec_t epsilon    = {1.0};
  int n_site       = 1;
  int n_bath       = epsilon.size();
  int n_spin       = cp.gf_struct.size();
  double U         = 0.0;
  double mu        = 0.0;
  double t         = 0.0;
  double tau_max   = cp.beta;
  double tau_split = cp.beta * 0.555;

  //parameters for TCI
  constexpr int n_GK = 15;
  auto [vi, wi]      = QuadratureGK<n_GK>(0, 1);
  int bondDim        = 10;
  int sweepBound     = 30;

  // prepare input
  auto [Delta_tau, ad_imp, u_tau, G_tau] = test_setup(n_site, n_bath, n_spin, U, mu, t, cp, theta, epsilon);
  auto u_interpolator                    = interpolator_t<scalar_t>(u_tau, u_tau[0].mesh().size());
  frame_t frame_zeroth_order             = u_interpolator(tau_max - tau_split) * u_interpolator(tau_split);
  long n_bl                              = cp.gf_struct.size();
  std::vector<int> block_shape           = {};
  for (auto [bl, bl_pair] : enumerate(cp.gf_struct)) {
    auto [bl_name, bl_size] = bl_pair;
    block_shape.push_back(bl_size);
  }
  std::vector<std::vector<fop_t>> all_d_ops     = {};
  std::vector<std::vector<fop_t>> all_d_dag_ops = {};
  all_d_ops.resize(n_bl, {});
  all_d_dag_ops.resize(n_bl, {});
  auto fops = fundamental_operator_set{cp.gf_struct};
  for (auto [bl, bl_pair] : enumerate(cp.gf_struct)) {
    auto [bl_name, bl_size] = bl_pair;
    all_d_ops[bl].clear();
    all_d_dag_ops[bl].clear();
    for (auto idx : range(bl_size)) {
      all_d_ops[bl].push_back({0.0, false, fops[{bl_name, idx}], bl, idx});
      all_d_dag_ops[bl].push_back({0.0, true, fops[{bl_name, idx}], bl, idx});
    }
  }
  std::cout << "Delta_tau shape:" << std::endl;
  print_block_shape(Delta_tau);
  std::cout << "G_tau shape:" << std::endl;
  print_block_shape(G_tau);
  std::cout << "u_tau shape:" << std::endl;
  print_block_shape(u_tau);
  std::cout << std::setprecision(17) << std::endl;

  // test decomposition of u_products_00
  bool do_naive               = false;
  bool do_split               = true;
  std::vector<int> order_list = {1, 2, 3, 4};
  std::vector<double> integral_list;
  if (do_naive) {
    for (int order : order_list) {
      int n = 2 * order;                      // number of tau's
      std::vector<int> iota_d_list(order, 0); // fix orbital indices temporarily
      std::vector<int> iota_d_dag_list(order, 0);
      double integral = 0.0;
      std::vector<int> pivot1(order, 0);
      pivot1.insert(pivot1.end(), order, n_GK - 2);

      std::vector<int> range(n);
      std::iota(range.begin(), range.end(), 0);
      auto order_list_pair  = get_all_order(range); //gives all possible phi
      long count            = 0;
      auto get_u_tau_max_00 = [&frame_zeroth_order, &tau_split, &tau_max, &all_d_ops, &all_d_dag_ops, &block_shape, &cp, &Delta_tau = Delta_tau,
                               &ad_imp = ad_imp, &u_interpolator, &count, &order_list_pair, &iota_d_list, &iota_d_dag_list](std::vector<double> vs) {
        auto config =
           BuildConfig(frame_zeroth_order, tau_split, tau_max, all_d_ops, all_d_dag_ops, block_shape, cp, Delta_tau, ad_imp, u_interpolator);
        auto taus    = changeVariable(vs, tau_max);
        double sum   = 0.0;
        double value = 0.0;
        for (auto [order_c_list, order_c_dag_list] : order_list_pair) {
          config(getElements(order_c_list, taus), getElements(order_c_dag_list, taus), iota_d_list, iota_d_dag_list);
          config.evaluate_hyb_weight();
          config.evaluate_u_products();
          if (config.u_products[0].size() == 0) {
            value = 0.0;
          } else {
            value = config.u_products[0](0, 0) * config.hyb_weight * config.sign;
          }
          sum += value;
        }
        count++;
        double j = jacobian(taus, tau_max);
        return sum * j;
      };
      std::vector<double> v1;
      for (int i = 0; i < pivot1.size(); i++) { v1.push_back(vi[pivot1[i]]); }
      std::cout << "v1: ";
      print_vector(v1);
      auto u_tau_max_00_v1 = get_u_tau_max_00(v1);
      std::cout << "get_u_tau_max_00(pivot1): " << u_tau_max_00_v1 << std::endl;

      /// do TCI
      // std::cout << "tci1" << std::endl;
      // auto ci1 = xfac::CTensorCI<double, double>(get_u_tau_max_00, std::vector(n, vi), {.pivot1 = pivot1});
      // for (int i = 0; i < sweepBound; i++) {
      //   ci1.iterate();
      //   integral = ci1.sumWeighted(std::vector(n, wi));
      //   std::cout << i << " " << count << " " << ci1.pivotError[ci1.pivotError.size() - 1] << " " << integral << std::endl;
      // }
      std::cout << "tci2" << std::endl;
      auto ci2 = xfac::CTensorCI2<double, double>(get_u_tau_max_00, std::vector(n, vi), {.bond_dim = bondDim, .pivot1 = pivot1});
      for (int i = 0; i < sweepBound; i++) {
        ci2.iterate();
        if (i == sweepBound - 1) { ci2.makeCanonical(); }
        integral = ci2.tt.sum(std::vector(n, wi));
        std::cout << i << " " << count << " " << ci2.pivotError[ci2.pivotError.size() - 1] << " " << integral << std::endl;
      }
      integral_list.push_back(integral);
    }
    std::cout << "u_tau_max_00 exact: " << u_interpolator(tau_max)[0](0, 0) << std::endl;
    std::cout << "order 0: " << frame_zeroth_order[0](0, 0) << std::endl;
    for (int i = 0; i < order_list.size(); i++) { std::cout << "order " << order_list[i] << ": " << integral_list[i] << std::endl; }
    std::cout << "TCI sum: " << frame_zeroth_order[0](0, 0) + std::accumulate(integral_list.begin(), integral_list.end(), 0.0) << std::endl;
  }

  if (do_split) {
    std::vector<double> calculation_time_list = {};
    for (int order : order_list) {
      auto start_time = std::chrono::high_resolution_clock::now();
      int n           = 2 * order;            // number of tau's
      std::vector<int> iota_d_list(order, 0); // fix orbital indices temporarily
      std::vector<int> iota_d_dag_list(order, 0);
      double integral = 0.0;
      for (int n_left = 1; n_left < n; n_left++) {
        int n_right = n - n_left;
        std::vector<int> pivot1_left(n_left, n_GK - 2);
        std::vector<int> pivot1_right(n_right, 0);
        std::vector<int> pivot1(pivot1_left);
        pivot1.insert(pivot1.end(), pivot1_right.begin(), pivot1_right.end());
        std::vector<int> range(n);
        std::iota(range.begin(), range.end(), 0);
        auto order_list_pair  = get_all_order(range); //gives all possible phi
        long count            = 0;
        auto get_u_tau_max_00 = [&frame_zeroth_order, &tau_split, &tau_max, &all_d_ops, &all_d_dag_ops, &block_shape, &cp, &Delta_tau = Delta_tau,
                                 &ad_imp = ad_imp, &u_interpolator, &count, &order_list_pair, &iota_d_list, &iota_d_dag_list,
                                 &n_left](const std::vector<double> &vs) {
          std::vector<double> vs_left(vs.begin(), vs.begin() + n_left);
          std::vector<double> vs_right(vs.begin() + n_left, vs.end());
          std::vector<double> taus_left  = changeVariable(vs_left, tau_split, 0.0);
          std::vector<double> taus_right = changeVariable(vs_right, tau_max, tau_split);
          auto taus(taus_left);
          taus.insert(taus.end(), taus_right.begin(), taus_right.end());
          // auto config =
          //    BuildConfig(frame_zeroth_order, tau_split, tau_max, all_d_ops, all_d_dag_ops, block_shape, cp, Delta_tau, ad_imp, u_interpolator);
          double sum   = 0.0;
          double value = 0.0;
          for (auto [order_c_list, order_c_dag_list] : order_list_pair) {
            // config(getElements(order_c_list, taus), getElements(order_c_dag_list, taus), iota_d_list, iota_d_dag_list);
            double value = evaluate_u_tau_max_00(frame_zeroth_order, tau_split, tau_max, all_d_ops, all_d_dag_ops, block_shape, cp, Delta_tau, ad_imp, u_interpolator, getElements(order_c_list, taus), getElements(order_c_dag_list, taus), iota_d_list, iota_d_dag_list);
            // config.evaluate_u_products();
            // if (config.u_products[0].size() == 0) {
            //   value = 0.0;
            // } else {
            //   value = config.u_products[0](0, 0) * config.hyb_weight * config.sign;
            // }
            sum += value;
          }
          count++;
          double j = jacobian(taus_left, tau_split, 0.0) * jacobian(taus_right, tau_max, tau_split);
          return sum * j;
        };
        std::vector<double> v1;
        for (int i = 0; i < pivot1.size(); i++) { v1.push_back(vi[pivot1[i]]); }
        std::cout << "pivot1: ";
        print_vector(v1);
        auto u_tau_max_00_v1 = get_u_tau_max_00(v1);
        std::cout << "get_u_tau_max_00(pivot1): " << u_tau_max_00_v1 << std::endl;
        if (u_tau_max_00_v1 == 0) { continue; }

        // do TCI
        // std::cout << "tci1" << std::endl;
        // auto ci1 = xfac::CTensorCI<double, double>(get_u_tau_max_00, std::vector(n, vi), {.pivot1 = pivot1});
        // double current_integral{0};
        // double previous_integral{0};
        // for (int i = 0; i < sweepBound; i++) {
        //   ci1.iterate();
        //   current_integral = ci1.sumWeighted(std::vector(n, wi));
        //   std::cout << i << " " << count << " " << ci1.pivotError[ci1.pivotError.size() - 1] << " " << current_integral << std::endl;
        //   if (std::abs(current_integral - previous_integral) < 1e-10 && i > 3) { break; }
        //   previous_integral = current_integral;
        // }
        std::cout << "tci2" << std::endl;
        auto ci2 = xfac::CTensorCI2<double, double>(get_u_tau_max_00, std::vector(n, vi), {.bond_dim = bondDim, .pivot1 = pivot1});
        double current_integral{0};
        double previous_integral{0};
        std::cout << "rank nEval LastSweepPivotError integral\n";
        for (int i = 0; i < sweepBound; i++) {
          ci2.iterate();
          if (i == sweepBound - 1) { ci2.makeCanonical(); }
          current_integral = ci2.tt.sum(std::vector(n, wi));
          std::cout << i << " " << count << " " << ci2.pivotError[ci2.pivotError.size() - 1] << " " << current_integral << std::endl;
          if (std::abs(current_integral - previous_integral) < 1e-8 && i > 1) { break; }
          previous_integral = current_integral;
        }
        std::cout << std::endl;
        integral += current_integral;
      }
      auto end_time            = std::chrono::high_resolution_clock::now();
      auto duration            = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
      auto duration_in_seconds = static_cast<double>(duration) / 1e6;
      calculation_time_list.push_back(duration_in_seconds);
      integral_list.push_back(integral);
    }
    std::cout << "u_tau_max_00 exact: " << std::setw(10) << u_interpolator(tau_max)[0](0, 0) << std::endl;

    std::cout << std::left << std::setw(10) << "order" << std::setw(30) << "value" << std::setw(30) << "time(s)" << std::endl;

    std::cout << std::setw(10) << "0" << std::setw(30) << frame_zeroth_order[0](0, 0) << std::endl;

    for (int i = 0; i < order_list.size(); i++) {
      std::cout << std::setw(10) << order_list[i] << std::setw(30) << integral_list[i] << std::setw(30) << calculation_time_list[i] << std::endl;
    }

    double sum_value = frame_zeroth_order[0](0, 0) + std::accumulate(integral_list.begin(), integral_list.end(), 0.0);
    double sum_time  = std::accumulate(calculation_time_list.begin(), calculation_time_list.end(), 0.0);

    std::cout << std::setw(10) << "sum:" << std::setw(30) << sum_value << std::setw(30) << sum_time << std::endl;

    return 0;
  }

  return 0;
}
