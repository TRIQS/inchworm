#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "./mode.hpp"

using namespace inchworm;

void ModeCombineFactorization::run_single_element() {

  // TCI
  int n_phi                                 = std::accumulate(mp.gf_block_shape.begin(), mp.gf_block_shape.end(), 0);
  for (int order : sp.order_list) {
    auto start_time = std::chrono::high_resolution_clock::now();
    int n           = 2 * order;     // number of tau's
    std::vector<int> v_pivot1(n, 0); // for tau only; pivots for iota are set later
    std::vector<int> range(n);
    std::iota(range.begin(), range.end(), 0);
    auto phi_pair_list      = get_all_phi(range);               //gives all possible phi
    auto iota_pair_list     = std::move(get_all_iota(mp.gf_block_shape, order)); //gives all possible iota
    double integral_sum_phi = 0.0;
    for (auto [phi_d_list, phi_d_dag_list] : phi_pair_list) {
      double integral_sum_n_left = 0.0;
      for (int n_left = 1; n_left < n; n_left++) {
        double integral_sum_iota = 0.0;
        // for (auto [iota_d_list, iota_d_dag_list] : iota_pair_list) {
        long count                 = 0;
        auto get_u_tau_max_element = [this,&count, &phi_d_list = phi_d_list,
                                      &phi_d_dag_list = phi_d_dag_list, &n_left, 
                                      &iota_pair_list](const std::vector<double> &v_iota_s) {
          std::vector<double> vs(v_iota_s.begin() + 1, v_iota_s.end());
          auto [iota_d_list, iota_d_dag_list] = iota_pair_list[static_cast<int>(v_iota_s[0])];
          std::vector<double> vs_left(vs.begin(), vs.begin() + n_left);
          std::vector<double> vs_right(vs.begin() + n_left, vs.end());
          std::vector<double> taus_left  = change_variable(vs_left, sp.tau_split, 0.0);
          std::vector<double> taus_right = change_variable(vs_right, sp.tau_max, sp.tau_split);
          auto taus(taus_left);
          taus.insert(taus.end(), taus_right.begin(), taus_right.end());
          double integrand = evaluate_u_tau_max(sr.u_tau_max_zeroth_order, sp.tau_split, sp.tau_max, mp.all_d_ops, mp.all_d_dag_ops, mp.gf_block_shape, cp, mp.Delta_tau,
                                                mp.ad_imp, sr.u_interpolator, get_elements(phi_d_list, taus), get_elements(phi_d_dag_list, taus), iota_d_list,
                                                iota_d_dag_list, sp.bl_index, sp.subspace_index);
          count++;
          double j = jacobian(taus_left, sp.tau_split, 0.0) * jacobian(taus_right, sp.tau_max, sp.tau_split);
          return integrand * j;
        };

        std::vector<double> v_iota_s1;
        for (int i = 0; i < v_pivot1.size(); i++) { v_iota_s1.push_back(tp.vi[v_pivot1[i]]); }
        double u_tau_max_element_vs1 = 0;
        // set pivot for iota
        std::vector<double> iota_pivot_list_valid = {};
        for (int iota_pivot1 = 0; iota_pivot1 < iota_pair_list.size(); iota_pivot1++) {
          auto v_iota_s1_temp = v_iota_s1;
          v_iota_s1_temp.insert(v_iota_s1_temp.begin(), iota_pivot1);
          u_tau_max_element_vs1 = get_u_tau_max_element(v_iota_s1_temp);
          if (u_tau_max_element_vs1 != 0) { iota_pivot_list_valid.push_back(iota_pivot1); }
        }
        if (iota_pivot_list_valid.size() == 0) {         std::cerr << "no valid iota found" << std::endl; continue; }
        else{ std::cout << "iota_pivot_list_valid size: " << iota_pivot_list_valid.size() << std::endl;
        std::cout << "iota_pair_list size: " << iota_pair_list.size() << std::endl;
        }
        v_iota_s1.insert(v_iota_s1.begin(), iota_pivot_list_valid[0]);

        auto pivot1 = v_pivot1;
        pivot1.insert(pivot1.begin(), 0);

        std::cout << "pivot1: ";
        print_vector(pivot1);
        std::cout << "v_iota_s1: ";
        print_vector(v_iota_s1);

        if (sp.debug) {
          std::vector<double> vs1(v_iota_s1.begin() + 1, v_iota_s1.end());
          auto [iota_d_list1, iota_d_dag_list1] = iota_pair_list[iota_pivot_list_valid[static_cast<int>(v_iota_s1[0])]];
          std::vector<double> vs1_left(vs1.begin(), vs1.begin() + n_left);
          std::vector<double> vs1_right(vs1.begin() + n_left, vs1.end());
          std::vector<double> taus1_left  = change_variable(vs1_left, sp.tau_split, 0.0);
          std::vector<double> taus1_right = change_variable(vs1_right, sp.tau_max, sp.tau_split);
          auto taus(taus1_left);
          taus.insert(taus.end(), taus1_right.begin(), taus1_right.end());
          std::cout << "iota_d_list: ";
          print_vector(iota_d_list1);
          std::cout << "iota_d_dag_list: ";
          print_vector(iota_d_dag_list1);
          std::cout << "tau_d_list: ";
          print_vector(get_elements(phi_d_list, taus));
          std::cout << "tau_d_dag_list: ";
          print_vector(get_elements(phi_d_dag_list, taus));
          auto u_tau_max_element_vs1 = get_u_tau_max_element(v_iota_s1);
          std::cout << "get_u_tau_max_element(pivot1): " << u_tau_max_element_vs1 << "\n" << std::endl;
        }
        if (u_tau_max_element_vs1 == 0) { continue; }

        double current_integral{0};
        double previous_integral{0};
        double integral_element{0};

        std::vector<double> iotai = iota_pivot_list_valid;
        auto wi_iota = std::vector(iota_pivot_list_valid.size(), 1.0);

        if (sp.debug) { std::cout << "iteration nEval LastSweepPivotError integral\n"; }
        auto input_to_append = std::vector(n, tp.vi);
        auto input           = std::vector(1, iotai);
        input.insert(input.end(), input_to_append.begin(), input_to_append.end());
        auto weight_to_append = std::vector(n,  tp.wi_v);
        auto weight           = std::vector(1, wi_iota);
        weight.insert(weight.end(), weight_to_append.begin(), weight_to_append.end());
        // std::cout << "input: " << std::endl;
        // for (auto v : input) { print_vector(v); }
        // std::cout << "weight: " << std::endl;
        // for (auto v : weight) { print_vector(v); }
        if (tp.tci_prrlu) {
          auto ci = xfac::CTensorCI2<double, double>(get_u_tau_max_element, input, {.bond_dim = tp.bond_dim, .pivot1 = pivot1});
          for (int i = 0; i <  tp.sweep_bound; i++) {
            ci.iterate();
            ci.makeCanonical();
            current_integral = ci.tt.sum(weight);
            if (sp.debug) { std::cout << i << " " << count << " " << ci.pivotError[ci.pivotError.size() - 1] << " " << current_integral << std::endl; }
            if (std::abs(current_integral - previous_integral) <  tp.error_bound && i > 1) { break; }
            previous_integral = current_integral;
          }
          integral_element = current_integral;
          if (sp.debug) { print_rank(ci.tt); }
        } else {
          auto ci = xfac::CTensorCI<double, double>(get_u_tau_max_element, input, {.pivot1 = pivot1});
          for (int i = 0; i <  tp.sweep_bound; i++) {
            ci.iterate();
            current_integral = ci.sumWeighted(weight);
            if (sp.debug) { std::cout << i << " " << count << " " << ci.pivotError[ci.pivotError.size() - 1] << " " << current_integral << std::endl; }
            if (std::abs(current_integral - previous_integral) <  tp.error_bound && i > 1) { break; }
            previous_integral = current_integral;
          }
          integral_element = current_integral;
          if (sp.debug) {
            std::cout << "rank:" << std::endl;
            print_vector(ci.rank());
          }
        }
        if (sp.debug) { std::cout << std::endl; }
        integral_sum_iota += integral_element;
        // }
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
