#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "./mode.hpp"

using namespace inchworm;

void ModeCombineFactorization::run_single_element() {


  for (int order : sp.order_list) {
    auto start_time = std::chrono::high_resolution_clock::now();
    int n           = 2 * order; //the number of tau's, i.e., the number of operators
    std::vector<int> v_pivot1(n, 0);
    std::vector<int> index_range(n);
    std::iota(index_range.begin(), index_range.end(), 0); //index_range is now 0,1,2 ... n-1
    auto phi_pair_list = get_all_phi(
       index_range); //gives all possible phi, which means after we generate taus and iotas, we need to use this to set the corresponding d or d^{\dagger}
    auto iota_pair_list     = get_all_iota(mp.gf_block_shape, order); //gives all possible iota
    double integral_sum_phi = 0.0;
    for (auto [phi_d_list, phi_d_dag_list] : phi_pair_list) {
      double integral_sum_n_left = 0.0;
      for (int n_left = 1; n_left < n; n_left++) {
        double integral_sum_iota   = 0.0;
        long count                 = 0;
        auto get_u_tau_max_element = [this, &count, &phi_d_list = phi_d_list, &phi_d_dag_list = phi_d_dag_list, &n_left,
                                      &iota_pair_list](const std::vector<double> &v_iota_s) {
          auto [iota_d_list, iota_d_dag_list] = iota_pair_list[static_cast<int>(v_iota_s[0])];
          std::vector<double> vs(v_iota_s.begin() + 1, v_iota_s.end());
          auto [taus_left,taus_right, taus] = obtain_taus(vs, n_left, sp.tau_split, sp.tau_max);
          double integrand = evaluate_u_tau_max(sr.u_tau_max_zeroth_order, sp.tau_split, sp.tau_max, mp.all_d_ops, mp.all_d_dag_ops,
                                                mp.gf_block_shape, cp, mp.Delta_tau, mp.ad_imp, sr.u_interpolator, get_elements(phi_d_list, taus),
                                                get_elements(phi_d_dag_list, taus), iota_d_list, iota_d_dag_list, sp.bl_index, sp.subspace_index);
          count++;
          double j = jacobian(taus_left, sp.tau_split, 0.0) * jacobian(taus_right, sp.tau_max, sp.tau_split);
          return integrand * j;
        };

        std::vector<double> v_iota_s1;
        for (int i = 0; i < v_pivot1.size(); i++) { v_iota_s1.push_back(tp.vi[v_pivot1[i]]); }
        double u_tau_max_element_vs1              = 0;
        std::vector<double> iota_pivot_list_valid = {}; // all the valid iota pivot
        for (int iota_pivot1 = 0; iota_pivot1 < iota_pair_list.size(); iota_pivot1++) {
          auto v_iota_s1_temp = v_iota_s1;
          v_iota_s1_temp.insert(v_iota_s1_temp.begin(), iota_pivot1);
          u_tau_max_element_vs1 = get_u_tau_max_element(v_iota_s1_temp);
          if (u_tau_max_element_vs1 != 0) { iota_pivot_list_valid.push_back(iota_pivot1); }
        }
        if (iota_pivot_list_valid.size() == 0) {
          if (sp.debug > 1) { std::cerr << "no valid iota found" << std::endl; }
          continue;
        } else {
          if (sp.debug > 1) {
            std::cout << "iota_pivot_list_valid size: " << iota_pivot_list_valid.size() << std::endl;
            std::cout << "iota_pair_list size: " << iota_pair_list.size() << std::endl;
          }
        }
        v_iota_s1.insert(v_iota_s1.begin(), iota_pivot_list_valid[0]);
        auto pivot1 = v_pivot1;
        pivot1.insert(pivot1.begin(), 0); // 0 refers to the first element of iota_pivot_list_valid

        u_tau_max_element_vs1 = get_u_tau_max_element(v_iota_s1);
        if (sp.debug > 1) {
          auto [iota_d_list1, iota_d_dag_list1] = iota_pair_list[iota_pivot_list_valid[static_cast<int>(v_iota_s1[0])]];
          std::vector<double> vs1(v_iota_s1.begin() + 1, v_iota_s1.end());
                   auto [taus_left1,taus_right1, taus1] = obtain_taus(vs1, n_left, sp.tau_split, sp.tau_max);
          print_pivot1(iota_d_list1, iota_d_dag_list1, get_elements(phi_d_list, taus1), get_elements(phi_d_dag_list, taus1), u_tau_max_element_vs1);
        }
        if (u_tau_max_element_vs1 == 0) { continue; }

        std::vector<double> iotai = iota_pivot_list_valid;
        auto wi_iota              = std::vector(iota_pivot_list_valid.size(), 1.0);
        auto input_to_append      = std::vector(n, tp.vi);
        auto input                = std::vector(1, iotai);
        input.insert(input.end(), input_to_append.begin(), input_to_append.end());
        auto weight_to_append = std::vector(n, tp.wi_v);
        auto weight           = std::vector(1, wi_iota);
        weight.insert(weight.end(), weight_to_append.begin(), weight_to_append.end());

        double integral_element = do_TCI<double, double>(get_u_tau_max_element, input, weight, pivot1, tp.sweep_bound, tp.bond_dim,
                                                         tp.integral_error_bound, tp.pivot_error_bound, tp.tci_prrlu, sp.debug, count);
        integral_sum_iota += integral_element;
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
