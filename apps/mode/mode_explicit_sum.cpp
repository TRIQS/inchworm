#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "./mode.hpp"

using namespace inchworm;

void ModeExplicitSum::run_single_element() {

  for (int order : sp.order_list) {
    auto start_time = std::chrono::high_resolution_clock::now();
    int n           = 2 * order; // number of tau's
    std::vector<int> pivot1(n, 0);
    std::vector<int> range(n);
    std::iota(range.begin(), range.end(), 0);
    auto phi_pair_list      = get_all_phi(range);               //gives all possible phi
    auto iota_pair_list     = get_all_iota(mp.gf_block_shape, order); //gives all possible iota
    double integral_sum_phi = 0.0;
    for (auto [phi_d_list, phi_d_dag_list] : phi_pair_list) {
      double integral_sum_n_left = 0.0;
      for (int n_left = 1; n_left < n; n_left++) {
        double integral_sum_iota = 0.0;
        for (auto [iota_d_list, iota_d_dag_list] : iota_pair_list) {
          long count                 = 0;
          auto get_u_tau_max_element = [this, &count, &phi_d_list = phi_d_list,
                                        &phi_d_dag_list = phi_d_dag_list, &iota_d_list = iota_d_list, &iota_d_dag_list = iota_d_dag_list, &n_left](const std::vector<double> &vs) {
            std::vector<double> vs_left(vs.begin(), vs.begin() + n_left);
            std::vector<double> vs_right(vs.begin() + n_left, vs.end());
            std::vector<double> taus_left  = change_variable(vs_left, sp.tau_split, 0.0);
            std::vector<double> taus_right = change_variable(vs_right, sp.tau_max, sp.tau_split);
            auto taus(taus_left);
            taus.insert(taus.end(), taus_right.begin(), taus_right.end());
            double integrand = evaluate_u_tau_max(sr.u_tau_max_zeroth_order, sp.tau_split, sp.tau_max, mp.all_d_ops, mp.all_d_dag_ops, mp.gf_block_shape, cp, mp.Delta_tau,
                                                  mp.ad_imp, sr.u_interpolator, get_elements(phi_d_list, taus), get_elements(phi_d_dag_list, taus),
                                                  iota_d_list, iota_d_dag_list, sp.bl_index, sp.subspace_index);
            count++;
            double j = jacobian(taus_left, sp.tau_split, 0.0) * jacobian(taus_right, sp.tau_max, sp.tau_split);
            return integrand * j;
          };

          std::vector<double> vs1;
          for (int i = 0; i < pivot1.size(); i++) { vs1.push_back(tp.vi[pivot1[i]]); }
          double u_tau_max_element_vs1 = 0;
          if (sp.debug) {
            std::vector<double> vs1_left(vs1.begin(), vs1.begin() + n_left);
            std::vector<double> vs1_right(vs1.begin() + n_left, vs1.end());
            std::vector<double> taus1_left  = change_variable(vs1_left, sp.tau_split, 0.0);
            std::vector<double> taus1_right = change_variable(vs1_right, sp.tau_max, sp.tau_split);
            auto taus(taus1_left);
            taus.insert(taus.end(), taus1_right.begin(), taus1_right.end());
            std::cout << "iota_d_list: ";
            print_vector(iota_d_list);
            std::cout << "iota_d_dag_list: ";
            print_vector(iota_d_dag_list);
            std::cout << "tau_d_list: ";
            print_vector(get_elements(phi_d_list, taus));
            std::cout << "tau_d_dag_list: ";
            print_vector(get_elements(phi_d_dag_list, taus));
            u_tau_max_element_vs1 = get_u_tau_max_element(vs1);
            std::cout << "get_u_tau_max_element(pivot1): " << u_tau_max_element_vs1 << "\n" << std::endl;
          } else {
            u_tau_max_element_vs1 = get_u_tau_max_element(vs1);
          }
          if (u_tau_max_element_vs1 == 0) { continue; }

          double current_integral{0};
          double previous_integral{0};
          double integral_element{0};
          if (sp.debug) { std::cout << "iteration nEval LastSweepPivotError integral\n"; }
          if (tp.tci_prrlu) {
            auto ci = xfac::CTensorCI2<double, double>(get_u_tau_max_element, std::vector(n, tp.vi), {.bond_dim = tp.bond_dim, .pivot1 = pivot1});
            for (int i = 0; i <  tp.sweep_bound; i++) {
              ci.iterate();
              ci.makeCanonical();
              current_integral = ci.tt.sum(std::vector(n,  tp.wi_v));
              if (sp.debug) { std::cout << i << " " << count << " " << ci.pivotError[ci.pivotError.size() - 1] << " " << current_integral << std::endl; }
              if (std::abs(current_integral - previous_integral) <  tp.error_bound && i > 1) { break; }
              previous_integral = current_integral;
            }
            integral_element = current_integral;
            if (sp.debug) {
              print_rank(ci.tt);
            }
          } else {
            auto ci = xfac::CTensorCI<double, double>(get_u_tau_max_element, std::vector(n, tp.vi), {.pivot1 = pivot1});
            for (int i = 0; i <  tp.sweep_bound; i++) {
              ci.iterate();
              current_integral = ci.sumWeighted(std::vector(n,  tp.wi_v));
              if (sp.debug) { std::cout << i << " " << count << " " << ci.pivotError[ci.pivotError.size() - 1] << " " << current_integral << std::endl; }
              if (std::abs(current_integral - previous_integral) <  tp.error_bound && i > 1) { break; }
              previous_integral = current_integral;
            }
            integral_element = current_integral;
          }
          if (sp.debug) { std::cout << std::endl; }
          integral_sum_iota += integral_element;
        }
        integral_sum_n_left += integral_sum_iota;
      }
      integral_sum_phi += integral_sum_n_left;
    }
    auto end_time            = std::chrono::high_resolution_clock::now();
    auto duration            = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    auto duration_in_seconds = static_cast<double>(duration) / 1e6;
    sr.calculation_time_list.push_back(duration_in_seconds);
    sr.integral_order_list.push_back(integral_sum_phi);
  }


}
