#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "./mode.hpp"

using namespace inchworm;

void ModePartitionFactorization::runSingleElement() {

  // TCI
  std::vector<double> integral_order_list   = {};
  std::vector<double> calculation_time_list = {};
  auto [vi, wi_v]                           = selectQuadratureGK(n_GK, 0, 1);
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
    auto phi_pair_list = getAllPhi(range); //gives all possible phi

    std::vector<int> iota_pivots(n, 0);
    std::vector<int> iota_pivots_range(n_phi);
    std::iota(iota_pivots_range.begin(), iota_pivots_range.end(), 0);
    std::vector<std::vector<int>> all_iota_pivots{};
    generateCombinations(iota_pivots_range, iota_pivots, 0, all_iota_pivots);
    // auto iota_pair_list     = getAllIota(block_shape, order); //gives all possible iota
    double integral_sum_phi = 0.0;
    for (auto [phi_d_list, phi_d_dag_list] : phi_pair_list) {
      double integral_sum_n_left = 0.0;
      for (int n_left = 1; n_left < n; n_left++) {
        double integral_sum_iota = 0.0;
        // for (auto [iota_d_list, iota_d_dag_list] : iota_pair_list) {
        long count                 = 0;
        auto get_u_tau_max_element = [this,  &count, &phi_d_list = phi_d_list,
                                      &phi_d_dag_list = phi_d_dag_list, &n_left](const std::vector<double> &v_iota_s) {
          int mid = v_iota_s.size() / 2;
          std::vector<double> iotas(v_iota_s.begin(), v_iota_s.begin() + mid);
          std::vector<double> vs(v_iota_s.begin() + mid, v_iota_s.end());
          int mid_iota = iotas.size() / 2;
          // std::cout << "vs: ";
          // printVector(vs);
          // std::cout << "iotas: ";
          // printVector(iotas);
          std::vector<double> iota_d_list(iotas.begin(), iotas.begin() + mid_iota);
          std::vector<double> iota_d_dag_list(iotas.begin() + mid_iota, iotas.end());
          std::vector<double> vs_left(vs.begin(), vs.begin() + n_left);
          std::vector<double> vs_right(vs.begin() + n_left, vs.end());
          std::vector<double> taus_left  = changeVariable(vs_left, tau_split, 0.0);
          std::vector<double> taus_right = changeVariable(vs_right, tau_max, tau_split);
          auto taus(taus_left);
          taus.insert(taus.end(), taus_right.begin(), taus_right.end());
          double integrand = evaluateUTauMax(u_tau_max_zeroth_order, tau_split, tau_max, all_d_ops, all_d_dag_ops, block_shape, cp, Delta_tau,
                                                ad_imp, u_interpolator, getElements(phi_d_list, taus), getElements(phi_d_dag_list, taus), iota_d_list,
                                                iota_d_dag_list, bl_index, subspace_index);
          count++;
          double j = jacobian(taus_left, tau_split, 0.0) * jacobian(taus_right, tau_max, tau_split);
          return integrand * j;
        };

        std::vector<double> v_iota_s1 {};
        double u_tau_max_element_vs1 = 0;

        // set pivot for iota
        int iota_pivot_index = 0;
        for (auto iota_pivot1 : all_iota_pivots) {
          std::vector<double> v_iota_s1_temp {};
          for (int i = 0; i < iota_pivot1.size(); i++) { v_iota_s1_temp.push_back(iotai[iota_pivot1[i]]); }
          for (int i = 0; i < v_pivot1.size(); i++) { v_iota_s1_temp.push_back(vi[v_pivot1[i]]); }
          u_tau_max_element_vs1 = get_u_tau_max_element(v_iota_s1_temp);
          if (u_tau_max_element_vs1 != 0) {
            v_iota_s1 = v_iota_s1_temp;
            break;
          }
          iota_pivot_index++;
        }
        if (iota_pivot_index == all_iota_pivots.size()) { continue; }

        auto pivot1           = all_iota_pivots[iota_pivot_index];
        auto pivot1_to_append = v_pivot1;
        pivot1.insert(pivot1.end(), pivot1_to_append.begin(), pivot1_to_append.end());

        std::cout << "pivot1: ";
        printVector(pivot1);
        std::cout << "v_iota_s1: ";
        printVector(v_iota_s1);

        if (debug) {
          int mid1 = v_iota_s1.size() / 2;
          std::vector<double> iotas1(v_iota_s1.begin(), v_iota_s1.begin() + mid1);
          std::vector<double> vs1(v_iota_s1.begin() + mid1, v_iota_s1.end());
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
          printVector(iota_d_list1);
          std::cout << "iota_d_dag_list: ";
          printVector(iota_d_dag_list1);
          std::cout << "tau_d_list: ";
          printVector(getElements(phi_d_list, taus));
          std::cout << "tau_d_dag_list: ";
          printVector(getElements(phi_d_dag_list, taus));
          std::cout << "get_u_tau_max_element(pivot1): " << u_tau_max_element_vs1 << "\n" << std::endl;
        }
        if (u_tau_max_element_vs1 == 0) { continue; }

        double current_integral{0};
        double previous_integral{0};
        double integral_element{0};
        if (debug) { std::cout << "iteration nEval LastSweepPivotError integral\n"; }
        auto input_to_append           = std::vector(n, vi);
        auto input = std::vector(n, iotai);
        input.insert(input.end(), input_to_append.begin(), input_to_append.end());
        auto weight_to_append           = std::vector(n, wi_v);
        auto weight = std::vector(n, wi_iota);
        weight.insert(weight.end(), weight_to_append.begin(), weight_to_append.end());
        // std::cout << "input: " << std::endl;
        // for (auto v : input) { printVector(v); }
        // std::cout << "weight: " << std::endl;
        // for (auto v : weight) { printVector(v); }
        if (tci_prrlu) {
          auto ci = xfac::CTensorCI2<double, double>(get_u_tau_max_element, input, {.bond_dim = bond_dim, .pivot1 = pivot1});
          for (int i = 0; i < sweep_bound; i++) {
            ci.iterate();
            ci.makeCanonical();
            current_integral = ci.tt.sum(weight);
            if (debug) { std::cout << i << " " << count << " " << ci.pivotError[ci.pivotError.size() - 1] << " " << current_integral << std::endl; }
            if (std::abs(current_integral - previous_integral) < error_bound && i > 1) { break; }
            previous_integral = current_integral;
          }
          integral_element = current_integral;
          if (debug) { printRank(ci.tt); }
        } else {
          auto ci = xfac::CTensorCI<double, double>(get_u_tau_max_element, input, {.pivot1 = pivot1});
          for (int i = 0; i < sweep_bound; i++) {
            ci.iterate();
            current_integral = ci.sumWeighted(weight);
            if (debug) { std::cout << i << " " << count << " " << ci.pivotError[ci.pivotError.size() - 1] << " " << current_integral << std::endl; }
            if (std::abs(current_integral - previous_integral) < error_bound && i > 1) { break; }
            previous_integral = current_integral;
          }
          integral_element = current_integral;
          if (debug) {
            std::cout << "rank:" << std::endl;
            printVector(ci.rank());
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
}
