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
    int n           = 2 * order; //the number of tau's, i.e., the number of operators
    std::vector<int> pivot1(n, 0);
    std::vector<int> index_range(n);
    std::iota(index_range.begin(), index_range.end(), 0); //index_range is now 0,1,2 ... n-1
    auto phi_pair_list      = get_all_phi(index_range);               //gives all possible phi, which means after we generate taus and iotas, we need to use this to set the corresponding d or d^{\dagger}
    auto iota_pair_list     = get_all_iota(mp.gf_block_shape, order); //gives all possible iota
    double integral_sum_phi = 0.0;
    for (auto [phi_d_list, phi_d_dag_list] : phi_pair_list) {
      double integral_sum_n_left = 0.0;
      for (int n_left = 1; n_left < n; n_left++) {
        double integral_sum_iota = 0.0;
        for (auto [iota_d_list, iota_d_dag_list] : iota_pair_list) {
          long count                 = 0;
          auto get_u_tau_max_element = [this, &count, &phi_d_list = phi_d_list,
                                        &phi_d_dag_list = phi_d_dag_list, &iota_d_list = iota_d_list, &iota_d_dag_list = iota_d_dag_list, &n_left](const std::vector<double> &vs)->double {
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
          double u_tau_max_element_vs1 = get_u_tau_max_element(vs1);
          if (sp.debug>1) {
            std::vector<double> vs1_left(vs1.begin(), vs1.begin() + n_left);
            std::vector<double> vs1_right(vs1.begin() + n_left, vs1.end());
            std::vector<double> taus1_left  = change_variable(vs1_left, sp.tau_split, 0.0);
            std::vector<double> taus1_right = change_variable(vs1_right, sp.tau_max, sp.tau_split);
            auto taus(taus1_left);
            taus.insert(taus.end(), taus1_right.begin(), taus1_right.end());
            print_pivot1(iota_d_list,iota_d_dag_list,get_elements(phi_d_list, taus),get_elements(phi_d_dag_list, taus),u_tau_max_element_vs1);
          } 
          if (u_tau_max_element_vs1 == 0) { continue; }
           

          auto input = std::vector(n, tp.vi);
          auto weight = std::vector(n, tp.wi_v);
          double integral_element = do_TCI<double,double> (get_u_tau_max_element, input, weight, pivot1, tp.sweep_bound, tp.bond_dim, tp.integral_error_bound, tp.pivot_error_bound, tp.tci_prrlu, sp.debug, count);
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
