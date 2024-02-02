#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "./mode.hpp"
#include <omp.h>

using namespace inchworm;

void ModeFullPartition::run_single_element() {

  for (int order : sp.order_list) {
    auto start_time = std::chrono::high_resolution_clock::now();
    int n           = 2 * order; //the number of tau's, i.e., the number of operators
    std::vector<int> pivot1(n, 0);
    std::vector<int> index_range(n);
    std::iota(index_range.begin(), index_range.end(), 0); //index_range is now 0,1,2 ... n-1
    // auto phi_pair_list = get_all_phi(
    //    index_range); //gives all possible phi, which means after we generate taus and iotas, we need to use this to set the corresponding d or d^{\dagger}
    auto iota_pair_list = get_all_iota(mp.gf_block_shape, order); //gives all possible iota
    std::vector<std::vector<int>> all_phi_pivots{};
    std::vector<int> phi_pivots_range{0, 1};
    std::vector<double> phi_value_range{0, 1};
    std::vector<int> phi_pivot(n, 0);
    generate_combinations(phi_pivots_range, phi_pivot, 0, all_phi_pivots);

    std::vector<std::vector<int>> all_n_left_pivots{};
    std::vector<int> n_left_pivots_range(n - 1);
    std::iota(n_left_pivots_range.begin(), n_left_pivots_range.end(), 0);
    std::vector<double> n_left_value_range(n - 1);
    std::iota(n_left_value_range.begin(), n_left_value_range.end(), 1);
    std::vector<int> n_left_pivot(1, 0);
    generate_combinations(n_left_pivots_range, n_left_pivot, 0, all_n_left_pivots);

    std::vector<std::vector<int>> all_iota_pivots{};
    std::vector<int> iota_pivots_range(mp.n_phi);
    std::iota(iota_pivots_range.begin(), iota_pivots_range.end(), 0);
    std::vector<double> iota_value_range(mp.n_phi);
    std::iota(iota_value_range.begin(), iota_value_range.end(), 0);
    std::vector<int> iota_pivot(n, 0);
    generate_combinations(iota_pivots_range, iota_pivot, 0, all_iota_pivots);

    double integral_sum_phi = 0.0;
    // for (auto [phi_d_list, phi_d_dag_list] : phi_pair_list) {
    double integral_sum_n_left = 0.0;
    double integral_sum_iota   = 0.0;
    // for (auto [iota_d_list, iota_d_dag_list] : iota_pair_list) {
    long count                 = 0;
    auto get_u_tau_max_element = [this, &count, &n](const std::vector<double> &fulls) -> double {
      int n_left = static_cast<int>(fulls[0]);
      // std::cout << "n_left = " << n_left << std::endl;
      std::vector<int> phi_list(fulls.begin() + 1, fulls.begin() + 1 + n);
      // std::cout << "phi_list = ";
      // print_vector(phi_list);
      std::vector<int> phi_d_list{};
      std::vector<int> phi_d_dag_list{};
      for (int i = 0; i < phi_list.size(); i++) {
        if (phi_list[i] == 0) {
          phi_d_list.push_back(i);
        } else if (phi_list[i] == 1) {
          phi_d_dag_list.push_back(i);
        } else {
          std::cout << "error in phi_list" << std::endl;
        }
      }
      // std::cout<< "phi_d_list = ";
      // print_vector(phi_d_list);
      // std::cout<< "phi_d_dag_list = ";
      // print_vector(phi_d_dag_list);
      if (phi_d_list.size() != phi_d_dag_list.size()) { return 0.0; }
      std::vector<int> iota_list(fulls.begin() + 1 + n, fulls.begin() + 1 + n + n);
      std::vector<int> iota_d_list           = get_elements(phi_d_list, iota_list);
      std::vector<int> iota_d_dag_list       = get_elements(phi_d_dag_list, iota_list);
      std::vector<int> number_in_block_d     = generate_number_in_block(mp.gf_block_shape, iota_d_list);
      std::vector<int> number_in_block_d_dag = generate_number_in_block(mp.gf_block_shape, iota_d_dag_list);
      if (number_in_block_d != number_in_block_d_dag) { return 0.0; }
      std::vector<double> vs(fulls.begin() + 1 + n + n, fulls.end());
      std::vector<double> vs_left(vs.begin(), vs.begin() + n_left);
      std::vector<double> vs_right(vs.begin() + n_left, vs.end());
      std::vector<double> taus_left  = change_variable(vs_left, sp.tau_split, 0.0);
      std::vector<double> taus_right = change_variable(vs_right, sp.tau_max, sp.tau_split);
      auto taus(taus_left);
      taus.insert(taus.end(), taus_right.begin(), taus_right.end());
      double integrand = evaluate_u_tau_max(sr.u_tau_max_zeroth_order, sp.tau_split, sp.tau_max, mp.all_d_ops, mp.all_d_dag_ops, mp.gf_block_shape,
                                            cp, mp.Delta_tau, mp.ad_imp, sr.u_interpolator, get_elements(phi_d_list, taus),
                                            get_elements(phi_d_dag_list, taus), iota_d_list, iota_d_dag_list, sp.bl_index, sp.subspace_index);
      count++;
      double j = jacobian(taus_left, sp.tau_split, 0.0) * jacobian(taus_right, sp.tau_max, sp.tau_split);
      return integrand * j;
    };

    std::vector<double> vs1;
    for (int i = 0; i < pivot1.size(); i++) { vs1.push_back(tp.v_value[pivot1[i]]); }
    std::vector<std::vector<int>> valid_pivots{};
    std::vector<int> valid_pivots_index    = {};
    std::vector<double> valid_pivots_value = {};
    bool found_valid_pivot                 = false;
// std::cout << "all_n_left_pivots.size() = " << all_n_left_pivots.size() << std::endl;
// std::cout << "all_phi_pivots.size() = " << all_phi_pivots.size() << std::endl;
// std::cout << "all_phi_pivots= " << std::endl;
// for (auto i : all_phi_pivots) {
//   print_vector(i);
//   std::cout << std::endl;
// }
#pragma omp parallel for collapse(3)
    for (auto pivot_n_left : all_n_left_pivots) {
      for (auto pivot_phi : all_phi_pivots) {
        for (auto pivot_iota : all_iota_pivots) {
          if (std::accumulate(pivot_phi.begin(), pivot_phi.end(), 0) == pivot_phi.size() / 2) {
            // std::cout << "pivot_n_left = ";
            // print_vector(pivot_n_left);
            // std::cout << "pivot_phi = ";
            // print_vector(pivot_phi);
            std::vector<double> fulls1;
            for (int i = 0; i < pivot_n_left.size(); i++) { fulls1.push_back(n_left_value_range[pivot_n_left[i]]); }
            for (int i = 0; i < pivot_phi.size(); i++) { fulls1.push_back(phi_value_range[pivot_phi[i]]); }
            for (int i = 0; i < pivot_iota.size(); i++) { fulls1.push_back(iota_value_range[pivot_iota[i]]); }
            fulls1.insert(fulls1.end(), vs1.begin(), vs1.end());
            // std::cout << "fulls1 = ";
            // print_vector(fulls1);
            // std::cout << "pivot1_temp = ";
            // print_vector(pivot1_temp);
            auto pivot_value = get_u_tau_max_element(fulls1);
            if (pivot_value != 0) {
              std::vector<int> pivot1_temp;
              for (int i = 0; i < pivot_n_left.size(); i++) { pivot1_temp.push_back(pivot_n_left[i]); }
              for (int i = 0; i < pivot_phi.size(); i++) { pivot1_temp.push_back(pivot_phi[i]); }
              for (int i = 0; i < pivot_iota.size(); i++) { pivot1_temp.push_back(pivot_iota[i]); }
              pivot1_temp.insert(pivot1_temp.end(), pivot1.begin(), pivot1.end());
#pragma omp critical
              {
                valid_pivots.push_back(pivot1_temp);
                valid_pivots_value.push_back(pivot_value);
                // std::cout << "pivot1_temp = ";
                // print_vector(pivot1_temp);
                // std::cout << "value: " << get_u_tau_max_element(fulls1) << std::endl;
                found_valid_pivot = true;
              }
            }
          }
        }
      }
    }
    std::cout << "size of valid_pivots = " << valid_pivots.size() << std::endl;
    // double u_tau_max_element_vs1 = get_u_tau_max_element(vs1);
    // if (sp.debug > 1) {
    //   auto [taus_left1, taus_right1, taus1] = obtain_taus(vs1, 1, sp.tau_split, sp.tau_max);
    //   print_pivot1(iota_d_list, iota_d_dag_list, get_elements(phi_d_list, taus1), get_elements(phi_d_dag_list, taus1), u_tau_max_element_vs1);
    // }
    if (!found_valid_pivot) { continue; }
    bool reduce_valid_pivots = true;
    if (reduce_valid_pivots) {
      valid_pivots_index.resize(valid_pivots.size());
      iota(valid_pivots_index.begin(), valid_pivots_index.end(), 0);
      sort_B_according_A(valid_pivots_value, valid_pivots_index);
      std::vector<std::vector<int>> valid_pivots_new{};
      for (int i = 0; i < valid_pivots_index.size(); i++) { valid_pivots_new.push_back(valid_pivots[valid_pivots_index[i]]); }
      valid_pivots = valid_pivots_new;
      std::cout << "size of valid_pivots (after reduction) = " << valid_pivots.size() << std::endl;
    }
    bool further_reduce = false;
    int max_size = 20;
    if(further_reduce){
      if(valid_pivots.size() > max_size){
        valid_pivots.resize(max_size);
      }
    }
    pivot1 = valid_pivots[0];

    auto input_n_left  = std::vector(1, n_left_value_range);
    auto weight_n_left = std::vector(1, std::vector<double>(n_left_value_range.size(), 1.0));
    auto input_phi     = std::vector(n, phi_value_range);
    auto weight_phi    = std::vector(n, std::vector<double>(phi_value_range.size(), 1.0));
    auto input_iota    = std::vector(n, iota_value_range);
    auto weight_iota   = std::vector(n, std::vector<double>(iota_value_range.size(), 1.0));
    auto input_v       = std::vector(n, tp.v_value);
    auto weight_v      = std::vector(n, tp.v_weight);
    auto input         = input_n_left;
    input.insert(input.end(), input_phi.begin(), input_phi.end());
    input.insert(input.end(), input_iota.begin(), input_iota.end());
    input.insert(input.end(), input_v.begin(), input_v.end());
    auto weight = weight_n_left;
    weight.insert(weight.end(), weight_phi.begin(), weight_phi.end());
    weight.insert(weight.end(), weight_iota.begin(), weight_iota.end());
    weight.insert(weight.end(), weight_v.begin(), weight_v.end());
    // std::cout << "input = ";
    // for(auto i:input){
    //   print_vector(i);
    //   std::cout << "\n";
    // }
    // std::cout << "weight = ";
    // for(auto i:weight){
    //   print_vector(i);
    //   std::cout << "\n";
    // }
    double integral_element =
       do_TCI_add_pivots<double, double>(get_u_tau_max_element, input, weight, pivot1, tp.sweep_bound, tp.bond_dim, tp.integral_error_bound,
                                         tp.pivot_error_bound, tp.tci_prrlu, sp.debug, count, valid_pivots);

    integral_sum_iota += integral_element;
    // }
    //for debugging
    integral_sum_n_left += integral_sum_iota;
    integral_sum_phi += integral_sum_n_left;
    // }
    auto end_time            = std::chrono::high_resolution_clock::now();
    auto duration            = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    auto duration_in_seconds = static_cast<double>(duration) / 1e6;
    sr.calculation_time_list.push_back(duration_in_seconds);
    sr.integral_list.push_back(integral_sum_phi);
  }
}
