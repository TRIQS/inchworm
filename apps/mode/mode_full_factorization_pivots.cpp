#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "./mode.hpp"

using namespace inchworm;

void ModeFullFactorizationPivots::run_single_element() {

  std::vector<double> iotai(mp.n_phi);
  std::iota(iotai.begin(), iotai.end(), 0);
  auto wi_iota = std::vector(mp.n_phi, 1.0);
  for (int order : sp.order_list) {
    auto start_time = std::chrono::high_resolution_clock::now();
    int n           = 2 * order;     // number of tau's
    std::vector<int> v_pivot1(n, 7); // for tau only; pivots for iota are set later
    std::vector<int> index_range(n);
    std::iota(index_range.begin(), index_range.end(), 0);
    auto phi_pair_list = get_all_phi(index_range); //gives all possible phi

    std::vector<int> iota_pivots(n, 0);           // this is an intermediate variable for generating all possible iota
    std::vector<int> iota_pivots_range(mp.n_phi); // the int version of iotai
    std::iota(iota_pivots_range.begin(), iota_pivots_range.end(), 0);
    std::vector<std::vector<int>> all_iota_pivots{};
    generate_combinations(iota_pivots_range, iota_pivots, 0, all_iota_pivots);

    double integral_sum_phi = 0.0;
    for (auto [phi_d_list, phi_d_dag_list] : phi_pair_list) {
      double integral_sum_n_left = 0.0;
      for (int n_left = 1; n_left < n; n_left++) {
        double integral_sum_iota   = 0.0;
        long count                 = 0;
        auto get_u_tau_max_element = [this, &count, &phi_d_list = phi_d_list, &phi_d_dag_list = phi_d_dag_list,
                                      &n_left](const std::vector<double> &v_iota_s) {
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
          std::vector<double> iota_d_list     = get_elements(phi_d_list, iotas);
          std::vector<double> iota_d_dag_list = get_elements(phi_d_dag_list, iotas);
          std::vector<int> iota_d_list_int(iota_d_list.begin(), iota_d_list.end());
          std::vector<int> iota_d_dag_list_int(iota_d_dag_list.begin(), iota_d_dag_list.end());
          std::vector<int> number_in_block_d     = generate_number_in_block(mp.gf_block_shape, iota_d_list_int);
          std::vector<int> number_in_block_d_dag = generate_number_in_block(mp.gf_block_shape, iota_d_dag_list_int);
          if (number_in_block_d != number_in_block_d_dag) { return 0.0; }
          auto [taus_left, taus_right, taus] = obtain_taus(vs, n_left, sp.tau_split, sp.tau_max);
          double integrand                   = evaluate_u_tau_max(sr.u_tau_max_zeroth_order, sp.tau_split, sp.tau_max, mp.all_d_ops, mp.all_d_dag_ops,
                                                                  mp.gf_block_shape, cp, mp.Delta_tau, mp.ad_imp, sr.u_interpolator, get_elements(phi_d_list, taus),
                                                                  get_elements(phi_d_dag_list, taus), iota_d_list, iota_d_dag_list, sp.bl_index, sp.subspace_index);
          count++;
          double j = jacobian(taus_left, sp.tau_split, 0.0) * jacobian(taus_right, sp.tau_max, sp.tau_split);
          return integrand * j;
        };

        std::vector<double> v_iota_s1{};
        std::vector<double> v_s1{};
        for (int i = 0; i < v_pivot1.size(); i++) { v_s1.push_back(tp.v_value[v_pivot1[i]]); } //insert v only temporarily
        double u_tau_max_element_vs1 = 0;

        // set pivot for iota
        int iota_pivot_index = 0;
        bool found_pivot1    = false;
        std::vector<int> valid_iota_index{};
        std::vector<double> valid_iota_value{};
        for (auto iota_pivot1 : all_iota_pivots) {
          std::vector<double> v_iota_s1_temp{};
          v_iota_s1_temp.reserve(v_s1.size() + iota_pivot1.size());
          for (int i = 0; i < iota_pivot1.size(); i++) {
            v_iota_s1_temp.push_back(iota_pivot1[i]);
            v_iota_s1_temp.push_back(v_s1[i]);
          }
          u_tau_max_element_vs1 = get_u_tau_max_element(v_iota_s1_temp);
          if (u_tau_max_element_vs1 != 0) {
            v_iota_s1    = v_iota_s1_temp;
            found_pivot1 = true;
            valid_iota_index.push_back(iota_pivot_index);
            valid_iota_value.push_back(u_tau_max_element_vs1);
          }
          iota_pivot_index++;
        }
        if (!found_pivot1) {
          // for debugging
          std::cout << "skipped" << std::endl;
          std::cout << "n_left: " << n_left << std::endl;
          std::cout << "phi_d_list: " << std::endl;
          for (auto i : phi_d_list) { std::cout << i << " "; }
          std::cout << std::endl;
          std::cout << "phi_d_dag_list: " << std::endl;
          for (auto i : phi_d_dag_list) { std::cout << i << " "; }
          std::cout << std::endl;
          continue;
        }
        std::cout << "valid_iota_index.size() before turncation: " << valid_iota_index.size() << std::endl;
        //truncate the valid_iota_index
        // std::cout << "valid_iota_value:";
        // print_vector(valid_iota_value);
        sort_B_according_A(valid_iota_value, valid_iota_index);
        if (std::abs(valid_iota_value[0]) < 1e-10) { continue; }
        //for debugging, print the size of valid_iota_index
        std::cout << "valid_iota_index.size(): " << valid_iota_index.size() << std::endl;
        std::vector<std::vector<int>> valid_pivots{};
        for (auto i : valid_iota_index) {
          std::vector<int> pivot1_temp = {};
          auto pivot1_to_append        = all_iota_pivots[i];
          pivot1_temp.reserve(v_pivot1.size() + pivot1_to_append.size());
          for (int i = 0; i < pivot1_to_append.size(); i++) {
            pivot1_temp.push_back(pivot1_to_append[i]);
            pivot1_temp.push_back(v_pivot1[i]);
          }
          valid_pivots.push_back(pivot1_temp);
        }
        auto pivot1 = valid_pivots[0];

        // std::vector<int> pivot1{};
        // auto pivot1_to_append = all_iota_pivots[iota_pivot_index];
        // pivot1.reserve(v_pivot1.size() + pivot1_to_append.size());
        // for (int i = 0; i < pivot1_to_append.size(); i++) {
        //   pivot1.push_back(pivot1_to_append[i]);
        //   pivot1.push_back(v_pivot1[i]);
        // }

        if (sp.debug > 1) {
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
          std::vector<double> iota_d_list1     = get_elements(phi_d_list, iotas1);
          std::vector<double> iota_d_dag_list1 = get_elements(phi_d_dag_list, iotas1);
          std::vector<int> iota_d_list_int1(iota_d_list1.begin(), iota_d_list1.end());
          std::vector<int> iota_d_dag_list_int1(iota_d_dag_list1.begin(), iota_d_dag_list1.end());
          auto [taus_left1, taus_right1, taus1] = obtain_taus(vs1, n_left, sp.tau_split, sp.tau_max);
          u_tau_max_element_vs1                 = valid_iota_value[0];
          print_pivot1(iota_d_list_int1, iota_d_dag_list_int1, get_elements(phi_d_list, taus1), get_elements(phi_d_dag_list, taus1),
                       u_tau_max_element_vs1);
        }
        if (u_tau_max_element_vs1 == 0) { continue; }

        std::vector<std::vector<double>> input{};
        input.reserve(2 * n);
        for (int i = 0; i < n; i++) {
          input.push_back(iotai);
          input.push_back(tp.v_value);
        }
        std::vector<std::vector<double>> weight{};
        weight.reserve(2 * n);
        for (int i = 0; i < n; i++) {
          weight.push_back(wi_iota);
          weight.push_back(tp.v_weight);
        }

        double integral_element =
           do_TCI_add_pivots<double, double>(get_u_tau_max_element, input, weight, pivot1, tp.sweep_bound, tp.bond_dim, tp.integral_error_bound,
                                             tp.pivot_error_bound, tp.tci_prrlu, sp.debug, count, valid_pivots);
        integral_sum_iota += integral_element;
        integral_sum_n_left += integral_sum_iota;
      }
      integral_sum_phi += integral_sum_n_left;
    }
    auto end_time            = std::chrono::high_resolution_clock::now();
    auto duration            = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    auto duration_in_seconds = static_cast<double>(duration) / 1e6;
    sr.calculation_time_list.push_back(duration_in_seconds);
    sr.integral_list.push_back(integral_sum_phi);
  }
}
