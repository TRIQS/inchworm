#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "./mode.hpp"

using namespace inchworm;

void ModeUseNormPivots::runSingleElement() {

  // TCI
  std::vector<double> integral_order_list   = {};
  std::vector<double> calculation_time_list = {};
  auto [vi, wi]                           = selectQuadratureGK(n_GK, 0, 1);
  for (int order : order_list) {
    auto start_time = std::chrono::high_resolution_clock::now();
    int n           = 2 * order; // number of tau's
    std::vector<int> pivot1(n, 0);
    std::vector<int> range(n);
    std::iota(range.begin(), range.end(), 0);
    auto phi_pair_list      = getAllPhi(range);               //gives all possible phi
    auto iota_pair_list     = getAllIota(block_shape, order); //gives all possible iota
    double integral_sum_phi = 0.0;
    for (auto [phi_d_list, phi_d_dag_list] : phi_pair_list) {
      double integral_sum_n_left = 0.0;
      for (int n_left = 1; n_left < n; n_left++) {
        double integral_sum_iota = 0.0;
        for (auto [iota_d_list, iota_d_dag_list] : iota_pair_list) {
          long count              = 0;
          auto get_u_tau_max_norm = [this, &count, &phi_d_list = phi_d_list, &phi_d_dag_list = phi_d_dag_list, &iota_d_list = iota_d_list,
                                     &iota_d_dag_list = iota_d_dag_list, &n_left](const std::vector<double> &vs) {
            std::vector<double> vs_left(vs.begin(), vs.begin() + n_left);
            std::vector<double> vs_right(vs.begin() + n_left, vs.end());
            std::vector<double> taus_left  = changeVariable(vs_left, tau_split, 0.0);
            std::vector<double> taus_right = changeVariable(vs_right, tau_max, tau_split);
            auto taus(taus_left);
            taus.insert(taus.end(), taus_right.begin(), taus_right.end());
            double integrand = evaluateUTauMax(u_tau_max_zeroth_order, tau_split, tau_max, all_d_ops, all_d_dag_ops, block_shape, cp, Delta_tau,
                                                  ad_imp, u_interpolator, getElements(phi_d_list, taus), getElements(phi_d_dag_list, taus),
                                                  iota_d_list, iota_d_dag_list, -1, subspace_index);
            count++;
            double j = jacobian(taus_left, tau_split, 0.0) * jacobian(taus_right, tau_max, tau_split);
            return integrand * j;
          };

          auto get_u_tau_max_element = [this, &count, &phi_d_list = phi_d_list, &phi_d_dag_list = phi_d_dag_list, &iota_d_list = iota_d_list,
                                        &iota_d_dag_list = iota_d_dag_list, &n_left](const std::vector<double> &vs) {
            std::vector<double> vs_left(vs.begin(), vs.begin() + n_left);
            std::vector<double> vs_right(vs.begin() + n_left, vs.end());
            std::vector<double> taus_left  = changeVariable(vs_left, tau_split, 0.0);
            std::vector<double> taus_right = changeVariable(vs_right, tau_max, tau_split);
            auto taus(taus_left);
            taus.insert(taus.end(), taus_right.begin(), taus_right.end());
            double integrand = evaluateUTauMax(u_tau_max_zeroth_order, tau_split, tau_max, all_d_ops, all_d_dag_ops, block_shape, cp, Delta_tau,
                                                  ad_imp, u_interpolator, getElements(phi_d_list, taus), getElements(phi_d_dag_list, taus),
                                                  iota_d_list, iota_d_dag_list, bl_index, subspace_index);
            count++;
            double j = jacobian(taus_left, tau_split, 0.0) * jacobian(taus_right, tau_max, tau_split);
            return integrand * j;
          };

          std::vector<double> vs1;
          for (int i = 0; i < pivot1.size(); i++) { vs1.push_back(vi[pivot1[i]]); }
          double u_tau_max_norm_vs1 = 0;
          if (debug) {
            std::vector<double> vs1_left(vs1.begin(), vs1.begin() + n_left);
            std::vector<double> vs1_right(vs1.begin() + n_left, vs1.end());
            std::vector<double> taus1_left  = changeVariable(vs1_left, tau_split, 0.0);
            std::vector<double> taus1_right = changeVariable(vs1_right, tau_max, tau_split);
            auto taus(taus1_left);
            taus.insert(taus.end(), taus1_right.begin(), taus1_right.end());
            std::cout << "iota_d_list: ";
            printVector(iota_d_list);
            std::cout << "iota_d_dag_list: ";
            printVector(iota_d_dag_list);
            std::cout << "tau_d_list: ";
            printVector(getElements(phi_d_list, taus));
            std::cout << "tau_d_dag_list: ";
            printVector(getElements(phi_d_dag_list, taus));
            u_tau_max_norm_vs1 = get_u_tau_max_norm(vs1);
            std::cout << "get_u_tau_max(pivot1): " << u_tau_max_norm_vs1 << std::endl;
            auto u_tau_max_element_vs1 = get_u_tau_max_element(vs1);
            std::cout << "get_u_tau_max_element(pivot1): " << u_tau_max_element_vs1 << "\n" << std::endl;
          } else {
            u_tau_max_norm_vs1 = get_u_tau_max_norm(vs1);
          }
          if (u_tau_max_norm_vs1 == 0) { continue; }

          double current_integral{0};
          double previous_integral{0};
          double integral_element{0};
          if (debug) { std::cout << "iteration nEval LastSweepPivotError integral\n"; }
          if (tci_prrlu) {
            auto ci = xfac::CTensorCI2<double, double>(get_u_tau_max_norm, std::vector(n, vi), {.bond_dim = bond_dim, .pivot1 = pivot1});
            for (int i = 0; i < sweep_bound; i++) {
              ci.iterate();
              ci.makeCanonical();
              current_integral = ci.tt.sum(std::vector(n, wi));
              if (debug) { std::cout << i << " " << count << " " << ci.pivotError[ci.pivotError.size() - 1] << " " << current_integral << std::endl; }
              if (std::abs(current_integral - previous_integral) < error_bound && i > 1) { break; }
              previous_integral = current_integral;
            }
            if (debug) { printRank(ci.tt); }
            auto u_tau_max_element_vs1 = get_u_tau_max_element(vs1);
            if (u_tau_max_element_vs1 != 0) {
              auto ci_element = xfac::CTensorCI2<double, double>(get_u_tau_max_element, std::vector(n, vi), {.pivot1 = pivot1});
              for (auto b = 0u; b < ci.len() - 1; b++) {
                auto pivots = ci.getPivotsAt(b);
                ci_element.addPivotsAt(pivots, b);
              }
              integral_element = ci_element.tt.sum(std::vector(n, wi));
            } else {
              integral_element = 0;
            }
          } else {
            auto ci = xfac::CTensorCI<double, double>(get_u_tau_max_norm, std::vector(n, vi), {.pivot1 = pivot1});
            for (int i = 0; i < sweep_bound; i++) {
              ci.iterate();
              current_integral = ci.sumWeighted(std::vector(n, wi));
              if (debug) { std::cout << i << " " << count << " " << ci.pivotError[ci.pivotError.size() - 1] << " " << current_integral << std::endl; }
              if (std::abs(current_integral - previous_integral) < error_bound && i > 1) { break; }
              previous_integral = current_integral;
            }
          }
          if (debug) { std::cout << std::endl; }
          integral_sum_iota += integral_element;
        }
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
