#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "./mode.hpp"

using namespace inchworm;

void ModePartitionBath::run_single_element() {

  std::vector<double> iotai(mp.n_phi);
  std::iota(iotai.begin(), iotai.end(), 0);
  auto wi_iota = std::vector(mp.n_phi, 1.0);
  for (int order : sp.order_list) {
    auto start_time = std::chrono::high_resolution_clock::now();
    int n           = 2 * order;     // number of tau's
    std::vector<int> v_pivot1(n, 0); // for tau only; pivots for iota are set later
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
        double integral_sum_iota        = 0.0;
        long count                      = 0;
        auto get_u_tau_max_element_init = [this, &count, &phi_d_list = phi_d_list, &phi_d_dag_list = phi_d_dag_list,
                                           &n_left](const std::vector<double> &v_iota_s) {
          int mid = v_iota_s.size() / 2;
          std::vector<double> iotas(v_iota_s.begin(), v_iota_s.begin() + mid);
          std::vector<double> vs(v_iota_s.begin() + mid, v_iota_s.end());
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
        double u_tau_max_element_vs1 = 0;

        // set pivot for iota
        int iota_pivot_index = 0;
        int index            = 0;
        bool found           = false;
        double max_value     = 0;
        for (auto iota_pivot1 : all_iota_pivots) {
          std::vector<double> v_iota_s1_temp{};
          for (int i = 0; i < iota_pivot1.size(); i++) { v_iota_s1_temp.push_back(iotai[iota_pivot1[i]]); }
          for (int i = 0; i < v_pivot1.size(); i++) { v_iota_s1_temp.push_back(tp.v_value[v_pivot1[i]]); }
          u_tau_max_element_vs1 = get_u_tau_max_element_init(v_iota_s1_temp);
          if (u_tau_max_element_vs1 != 0 && found == false) {
            iota_pivot_index = index;
            v_iota_s1        = v_iota_s1_temp;
            found            = true;
          }
          if (std::abs(u_tau_max_element_vs1) > std::abs(max_value)) { max_value = u_tau_max_element_vs1; }
          index++;
        }
        if (!found) {
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

        double bath_value               = max_value / 2;
        auto pivot1           = all_iota_pivots[iota_pivot_index];
        auto pivot1_to_append = v_pivot1;
        pivot1.insert(pivot1.end(), pivot1_to_append.begin(), pivot1_to_append.end());

        u_tau_max_element_vs1 = get_u_tau_max_element_init(v_iota_s1);
        if (sp.debug > 1) {
          int mid1 = v_iota_s1.size() / 2;
          std::vector<double> iotas1(v_iota_s1.begin(), v_iota_s1.begin() + mid1);
          std::vector<double> vs1(v_iota_s1.begin() + mid1, v_iota_s1.end());
          std::vector<double> iota_d_list1     = get_elements(phi_d_list, iotas1);
          std::vector<double> iota_d_dag_list1 = get_elements(phi_d_dag_list, iotas1);
          std::vector<int> iota_d_list_int1(iota_d_list1.begin(), iota_d_list1.end());
          std::vector<int> iota_d_dag_list_int1(iota_d_dag_list1.begin(), iota_d_dag_list1.end());
          auto [taus_left1, taus_right1, taus1] = obtain_taus(vs1, n_left, sp.tau_split, sp.tau_max);
          print_pivot1(iota_d_list_int1, iota_d_dag_list_int1, get_elements(phi_d_list, taus1), get_elements(phi_d_dag_list, taus1),
                       u_tau_max_element_vs1);
          std::cout << "bath value: " << bath_value << std::endl;
        }
        if (u_tau_max_element_vs1 == 0) { continue; }
        auto get_u_tau_max_element = [this, &count, &phi_d_list = phi_d_list, &phi_d_dag_list = phi_d_dag_list,
                                           &n_left,&bath_value](const std::vector<double> &v_iota_s) {
          for (int i = 0; i < v_iota_s.size(); i++) {
            if (v_iota_s[i] < 0) { return bath_value; }
          }
          int mid = v_iota_s.size() / 2;
          std::vector<double> iotas(v_iota_s.begin(), v_iota_s.begin() + mid);
          std::vector<double> vs(v_iota_s.begin() + mid, v_iota_s.end());
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
        

        auto vi_bath = tp.v_value;
        vi_bath.push_back(-1);
        auto iotai_bath = iotai;
        iotai_bath.push_back(-1);
        auto wi_bath = tp.v_weight;
        wi_bath.push_back(0.0);
        auto wi_iota_bath = wi_iota;
        wi_iota_bath.push_back(0.0);
        auto input_to_append = std::vector(n, vi_bath);
        auto input           = std::vector(n, iotai_bath);
        input.insert(input.end(), input_to_append.begin(), input_to_append.end());
        auto weight_to_append = std::vector(n, wi_bath);
        auto weight           = std::vector(n, wi_iota_bath);
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
    sr.integral_list.push_back(integral_sum_phi);
  }
}
