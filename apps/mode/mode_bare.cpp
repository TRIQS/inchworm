#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "./mode.hpp"

using namespace inchworm;

void ModeBare::run_single_element() {}

void ModeBare::run() {
  set_up();
  evaluate_partition_function();
  evaluate_greens_function();
}

void ModeBare::set_up() {}

void ModeBare::evaluate_partition_function() {
  // set up the iota values and weights
  std::vector<double> iota_value(mp.n_phi);
  std::iota(iota_value.begin(), iota_value.end(),
            0); // the value of iota is 0, 1, 2, 3, ... n_phi-1 where n_phi is the number of available spin-orbital labels
  std::vector<double> iota_weight(mp.n_phi, 1.0);
  for (int order : sp.order_list) {
    // initialization
    std::cout << "### order " << order << " ###" << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    int n           = 2 * order; // number of operators

    // prepare input, weight, and initial pivot for pretraining
    auto tci_input_pre        = std::vector(n, iota_value);
    auto tci_input_pre_append = std::vector(n, std::vector<double>{tp.v_value[int(tp.v_value.size() / 2)]});
    tci_input_pre.insert(tci_input_pre.end(), tci_input_pre_append.begin(), tci_input_pre_append.end());
    auto tci_weight_pre        = std::vector(n, iota_weight);
    auto tci_weight_pre_append = std::vector(n, std::vector<double>{1.0});
    tci_weight_pre.insert(tci_weight_pre.end(), tci_weight_pre_append.begin(), tci_weight_pre_append.end());
    std::vector<int> tci_pivot_pre{};
    tci_pivot_pre.insert(tci_pivot_pre.end(), n, 1); // the first part of the pivot is set to 1 to avoid 0 for the pretraining step
    tci_pivot_pre.insert(tci_pivot_pre.end(), n, 0);

    // prepare input, weight for training
    // the initial pivot for training is determined by the pretraining step
    auto tci_input        = std::vector(n, iota_value);
    auto tci_input_append = std::vector(n, tp.v_value);
    tci_input.insert(tci_input.end(), tci_input_append.begin(), tci_input_append.end());
    auto tci_weight        = std::vector(n, iota_weight);
    auto tci_weight_append = std::vector(n, tp.v_weight);
    tci_weight.insert(tci_weight.end(), tci_weight_append.begin(), tci_weight_append.end());
    std::vector<int> tci_pivot{};

    // generate all possible ordering of creation/annihiation operators
    // essentially, this is picking up n/2 operators from n operators
    std::vector<int> index_range(n);
    std::iota(index_range.begin(), index_range.end(), 0);
    auto phi_pair_list = get_all_phi(
       index_range); // this gives all binom(n, n/2) possible combinations; each element is a pair of vectors indicating the time ordering, one for creation and one for annihilation
    int n_phi_pair = phi_pair_list.size();

    // calculate the lower bound of the pre-trained integral
    // tp.integral_lower_bound is the bound for the integral contribution for a specific order, below which the integral can be skipped
    auto max_weight_v               = *std::max_element(tp.v_weight.begin(), tp.v_weight.end());
    double pre_integral_lower_bound = tp.integral_lower_bound / (std::pow(max_weight_v * tp.v_weight.size(), n) * n_phi_pair * (n - 1));
    if (sp.debug > 1) { std::cout << "pre_integral_lower_bound: " << pre_integral_lower_bound << std::endl; }

    // calculate the integral for the auxiliary function for the pre-training step
    int count_auxiliary = 0;
    auto get_auxiliary  = [this, &count_auxiliary](const std::vector<double> &v_iota_s) {
      int mid = v_iota_s.size() / 2;
      std::vector<double> iotas(v_iota_s.begin(), v_iota_s.begin() + mid);
      double sin_term = sin_func_all({}, iotas, mp.n_phi);
      count_auxiliary++;
      return tp.auxi_height * sin_term;
    };
    auto ci_pre_auxiliary = xfac::CTensorCI2<double, double>(
       get_auxiliary, tci_input_pre, {.bondDim = tp.bond_dim, .reltol = tp.reltol, .pivot1 = tci_pivot_pre, .fullPiv = false});
    int ci_count_auxiliary                = 0;
    double previous_pivot_error_auxiliary = -1E5;
    if (sp.debug > 1) {
      std::cout << "integral for the auxiliary function" << std::endl;
      std::cout << "iteration nEval LastSweepPivotError\n";
    }
    while (true) {
      ci_pre_auxiliary.iterate();
      auto last_pivot_error = ci_pre_auxiliary.pivotError[ci_pre_auxiliary.pivotError.size() - 1];
      ci_count_auxiliary++;
      if (std::abs(last_pivot_error - previous_pivot_error_auxiliary) < tp.convergence_bound && ci_count_auxiliary > tp.convergence_iter) { break; }
      previous_pivot_error_auxiliary = last_pivot_error;
      if (sp.debug > 1) {
        std::cout << ci_count_auxiliary << " " << ci_pre_auxiliary.tt.sum(tci_weight_pre) << " " << last_pivot_error << std::endl;
        print_rank(ci_pre_auxiliary.tt);
      }
    }
    double integral_auxiliary = ci_pre_auxiliary.tt.sum(tci_weight_pre);

    //step up finished
    //starting pretraining and training for each phi
    double integral_sum_phi      = 0.0;
    double time_for_pre_training = 0.0;
    double time_for_find_pivot   = 0.0;
    double time_for_training     = 0.0;
    for (auto [phi_d_list, phi_d_dag_list] : phi_pair_list) {
      if (sp.debug > 1) { std::cout << " ------- tci start ------- " << std::endl; }
      // define function for pretraining and pretraining
      long count_pre             = 0;
      auto integrand_pretraining = [this, &count_pre, &phi_d_list = phi_d_list,
                                    &phi_d_dag_list = phi_d_dag_list](const std::vector<double> &v_iota_s) {
        int mid = v_iota_s.size() / 2;
        std::vector<double> iotas(v_iota_s.begin(), v_iota_s.begin() + mid);
        std::vector<double> vs(v_iota_s.begin() + mid, v_iota_s.end());
        std::vector<double> iota_d_list     = get_elements(phi_d_list, iotas);
        std::vector<double> iota_d_dag_list = get_elements(phi_d_dag_list, iotas);
        std::vector<int> iota_d_list_int(iota_d_list.begin(), iota_d_list.end());
        std::vector<int> iota_d_dag_list_int(iota_d_dag_list.begin(), iota_d_dag_list.end());
        std::vector<int> number_in_block_d     = generate_number_in_block(mp.gf_block_shape, iota_d_list_int);
        std::vector<int> number_in_block_d_dag = generate_number_in_block(mp.gf_block_shape, iota_d_dag_list_int);
        double sin_term                        = sin_func_all({}, iotas, mp.n_phi);

        if (number_in_block_d != number_in_block_d_dag) { return tp.auxi_height * sin_term; }

        auto [taus_left, taus_right, taus] = obtain_taus(vs, 0, 0.0, cp.beta);

        double integrand = evaluate_u_tau_max(sr.u_tau_max_zeroth_order_bare, 0.0, cp.beta, mp.all_d_ops, mp.all_d_dag_ops, mp.gf_block_shape, cp,
                                              mp.Delta_tau, mp.ad_imp, {}, get_elements(phi_d_list, taus), get_elements(phi_d_dag_list, taus),
                                              iota_d_list, iota_d_dag_list, -1, -1);

        count_pre++;
        double j = jacobian(taus, cp.beta, 0.0);
        // std::cout << "integrand: " << integrand << std::endl;
        return integrand * j + tp.auxi_height * sin_term;
      };
      long count              = 0;
      auto integrand_training = [this, &count, &phi_d_list = phi_d_list, &phi_d_dag_list = phi_d_dag_list](const std::vector<double> &v_iota_s) {
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
        auto [taus_left, taus_right, taus] = obtain_taus(vs, 0, 0.0, cp.beta);
        double integrand = evaluate_u_tau_max(sr.u_tau_max_zeroth_order_bare, 0.0, cp.beta, mp.all_d_ops, mp.all_d_dag_ops, mp.gf_block_shape, cp,
                                              mp.Delta_tau, mp.ad_imp, {}, get_elements(phi_d_list, taus), get_elements(phi_d_dag_list, taus),
                                              iota_d_list, iota_d_dag_list, -1, -1);
        count++;
        double j = jacobian(taus, cp.beta, 0.0);
        return integrand * j;
      };

      // pretraining
      auto start_time_pre_training    = std::chrono::high_resolution_clock::now();
      auto ci_pre                     = xfac::CTensorCI2<double, double>(integrand_pretraining, tci_input_pre,
                                                     {.bondDim = tp.bond_dim, .reltol = tp.reltol, .pivot1 = tci_pivot_pre, .fullPiv = false});
      int ci_count_pre                = 0;
      double previous_pivot_error_pre = -1E5;
      if (sp.debug > 1) {
        std::cout << "integral for the pretraining function" << std::endl;
        std::cout << "iteration nEval LastSweepPivotError\n";
      }
      while (true) {
        ci_pre.iterate();
        auto last_pivot_error = ci_pre.pivotError[ci_pre.pivotError.size() - 1];
        ci_count_pre++;
        if (std::abs(last_pivot_error - previous_pivot_error_pre) < tp.convergence_bound && ci_count_pre > tp.convergence_iter) { break; }
        previous_pivot_error_pre = last_pivot_error;
        if (sp.debug > 1) {
          std::cout << ci_count_pre << " " << ci_pre.tt.sum(tci_weight_pre) << " " << last_pivot_error << std::endl;
          print_rank(ci_pre.tt);
        }
      }
      double integral_pre        = ci_pre.tt.sum(tci_weight_pre);
      auto end_time_pre_training = std::chrono::high_resolution_clock::now();
      auto duration_pre_training = std::chrono::duration_cast<std::chrono::microseconds>(end_time_pre_training - start_time_pre_training).count();
      auto duration_in_seconds_pre_training = static_cast<double>(duration_pre_training) / 1e6;

      double integral_diff = integral_pre - integral_auxiliary;
      if (sp.debug > 1) {
        std::cout << "duration_pre_training: " << duration_in_seconds_pre_training << std::endl;
        std::cout << "integral_pre: " << integral_pre << std::endl;
        std::cout << "integral_auxiliary: " << integral_auxiliary << std::endl;
        std::cout << "integral_pre-integral_auxiliary: " << integral_diff << std::endl;
      }
      time_for_pre_training += duration_in_seconds_pre_training;
      if (std::abs(integral_diff) < pre_integral_lower_bound) {
        std::cout << "pre_trained integral is too small, skip the integral" << std::endl;
        std::cout << "integral_diff: " << integral_diff << std::endl;
        std::cout << " ------- tci finish------- " << std::endl;
        continue;
      }

      auto start_time_find_pivot = std::chrono::high_resolution_clock::now();
      if (sp.debug > 1) {
        std::cout << "integral for the training function" << std::endl;
        std::cout << "iteration nEval LastSweepPivotError\n";
      }
      bool find_pivot1              = false;
      std::vector<int> pivot1_train = {};
      for (auto b = 0u; b < ci_pre.len() - 1; b++) {
        auto pivots = ci_pre.getPivotsAt(b);
        for (auto p : pivots) {
          auto mid = p.size() / 2;
          std::vector<double> v_iota_s(p.begin(), p.end() - mid);
          for (int i = 0; i < n; i++) { v_iota_s.push_back(tp.v_value[int(tp.v_value.size() / 2)]); }
          auto val = integrand_training(v_iota_s);
          if (val != 0) {
            find_pivot1  = true;
            pivot1_train = p;
            if (sp.debug > 1) { std::cout << "found pivot1 with value: " << val << std::endl; }
            break;
          }
        }
        if (find_pivot1) { break; }
      }
      if (!find_pivot1) {
        std::cout << "weird!! " << std::endl;
        std::cout << " ------- tci finish ------- " << std::endl;
        continue;
      }
      auto end_time_find_pivot = std::chrono::high_resolution_clock::now();
      auto duration_find_pivot = std::chrono::duration_cast<std::chrono::microseconds>(end_time_find_pivot - start_time_find_pivot).count();
      auto duration_in_seconds_find_pivot = static_cast<double>(duration_find_pivot) / 1e6;
      time_for_find_pivot += duration_in_seconds_find_pivot;
      if (sp.debug > 1) { std::cout << "iteration nEval LastSweepPivotError integral\n"; }
      auto start_time_training = std::chrono::high_resolution_clock::now();
      auto ci                  = xfac::CTensorCI2<double, double>(integrand_training, tci_input,
                                                 {.bondDim = tp.bond_dim, .reltol = tp.reltol, .pivot1 = pivot1_train, .fullPiv = false});
      ci.addPivots(ci_pre);
      if (sp.debug > 1) {
        print_rank(ci.tt);
        std::cout << "bond_dim: " << ci.param.bondDim << std::endl;
      }
      double current_integral = 0.0;
      for (int i = 1; i <= tp.sweep_bound; i++) {
        ci.iterate();
        current_integral      = ci.tt.sum(tci_weight);
        auto last_pivot_error = ci.pivotError[ci.pivotError.size() - 1];
        std::cout << i << " " << count << " " << last_pivot_error << " " << current_integral << std::endl;
        print_rank(ci.tt);
      }
      auto integral_element             = current_integral;
      auto end_time_training            = std::chrono::high_resolution_clock::now();
      auto duration_training            = std::chrono::duration_cast<std::chrono::microseconds>(end_time_training - start_time_training).count();
      auto duration_in_seconds_training = static_cast<double>(duration_training) / 1e6;
      time_for_training += duration_in_seconds_training;
      integral_sum_phi += integral_element;
      if (sp.debug > 1) { std::cout << "integral_element: " << integral_element << std::endl; }
      std::cout << " ------- tci finish ------- " << std::endl;
    } // end of phi loop
    auto end_time            = std::chrono::high_resolution_clock::now();
    auto duration            = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    auto duration_in_seconds = static_cast<double>(duration) / 1e6;
    sr.calculation_time_list.push_back(duration_in_seconds);
    sr.integral_list.push_back(integral_sum_phi);
    sr.pretrain_time_list.push_back(time_for_pre_training);
    sr.find_pivot_time_list.push_back(time_for_find_pivot);
    sr.train_time_list.push_back(time_for_training);
  } // end of order loop
}

void ModeBare::evaluate_greens_function() {}