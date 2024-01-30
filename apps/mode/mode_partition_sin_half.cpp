#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "./mode.hpp"

using namespace inchworm;

void ModePartitionSinHalf::run_single_element() {

  std::vector<double> iotai(mp.n_phi);
  std::iota(iotai.begin(), iotai.end(), 0);
  auto wi_iota                        = std::vector(mp.n_phi, 1.0);
  auto max_weight_v                   = *std::max_element(tp.wi_v.begin(), tp.wi_v.end());
  double auxi_height                  = tp.auxi_height; //height of the auxiliary function in the pre-training
  double reltol_test                  = tp.reltol;
  double tci_convergence_threshold    = 1E-20;
  double pre_train_relative_threshold = 1E-20;
  double time_for_find_pivot          = 0.0;
  double time_for_pre_training        = 0.0;
  for (int order : sp.order_list) {
    auto start_time = std::chrono::high_resolution_clock::now();
    int n           = 2 * order; // number of operators
    std::cout << "### order " << order << " ###" << std::endl;
    std::vector<int> v_pivot1(n, int(tp.vi.size() / 2)); // for tau only; pivots for iota are set later
    std::vector<int> v_pivot1_pre(n, 0);

    std::vector<int> index_range(n);
    std::iota(index_range.begin(), index_range.end(), 0);
    auto phi_pair_list = get_all_phi(index_range); //gives all possible phi
    int n_phi_pair     = phi_pair_list.size();

    std::vector<int> iota_pivots(n, 0);           // this is an intermediate variable for generating all possible iota
    std::vector<int> iota_pivots_range(mp.n_phi); // the int version of iotai
    std::iota(iota_pivots_range.begin(), iota_pivots_range.end(), 0);
    std::vector<std::vector<int>> all_iota_pivots{};
    generate_combinations(iota_pivots_range, iota_pivots, 0, all_iota_pivots);

    // calculate the lower bound of the pre-trained integral
    // tp.integral_lower_bound is the bound for the integral contribution for a specific order, below which the integral can be skipped
    double pre_integral_lower_bound = tp.integral_lower_bound / (std::pow(max_weight_v * tp.wi_v.size(), n) * n_phi_pair * (n - 1));
    // this bound is underestimated, ajust for now
    pre_integral_lower_bound = pre_integral_lower_bound * 10;
    std::cout << "pre_integral_lower_bound: " << pre_integral_lower_bound << std::endl;

    // calculate the integral for the auxiliary function for the pre-training
    int count_auxiliary = 0;
    auto get_auxiliary  = [this, &auxi_height, &count_auxiliary](const std::vector<double> &v_iota_s) {
      int mid = v_iota_s.size() / 2;
      std::vector<double> iotas(v_iota_s.begin(), v_iota_s.begin() + mid);
      std::vector<double> vs(v_iota_s.begin() + mid, v_iota_s.end());
      double sin_term = sin_func_all({}, iotas, mp.n_phi);
      // double sin_term = linear_func_all( iotas, mp.n_phi);
      count_auxiliary++;
      return auxi_height * sin_term;
    };
    auto pivot1_auxiliary = std::vector<int>(n, 1);
    pivot1_auxiliary.insert(pivot1_auxiliary.end(), v_pivot1_pre.begin(), v_pivot1_pre.end());
    auto input_to_append_pre  = std::vector(n, std::vector<double>{tp.vi[int(tp.vi.size() / 2)]}); // for the pre-training, the v variable is fixed
    auto input_pre            = std::vector(n, iotai);
    auto weight_artificial    = std::vector<double>{1.0};
    auto weight_to_append_pre = std::vector(n, weight_artificial);
    auto weight_pre           = std::vector(n, wi_iota);
    input_pre.insert(input_pre.end(), input_to_append_pre.begin(), input_to_append_pre.end());
    weight_pre.insert(weight_pre.end(), weight_to_append_pre.begin(), weight_to_append_pre.end());

    auto ci_pre_auxiliary = xfac::CTensorCI2<double, double>(
       get_auxiliary, input_pre, {.bondDim = tp.bond_dim, .reltol = reltol_test, .pivot1 = pivot1_auxiliary, .fullPiv = true});
    std::cout << "integral for the auxiliary function" << std::endl;
    std::cout << "iteration nEval LastSweepPivotError\n";
    int ci_count_auxiliary                = 0;
    double previous_pivot_error_auxiliary = -1E5;
    while (true) {
      ci_pre_auxiliary.iterate();
      auto last_pivot_error = ci_pre_auxiliary.pivotError[ci_pre_auxiliary.pivotError.size() - 1];
      std::cout << ci_count_auxiliary << " " << count_auxiliary << " " << last_pivot_error << " " << std::endl;
      print_rank(ci_pre_auxiliary.tt);
      ci_count_auxiliary++;
      if (std::abs(last_pivot_error - previous_pivot_error_auxiliary) < tci_convergence_threshold && ci_count_auxiliary > 3) { break; }
      previous_pivot_error_auxiliary = last_pivot_error;
    }
    double integral_auxiliary = ci_pre_auxiliary.tt.sum(weight_pre);
    std::cout << "integral_auxiliary: " << integral_auxiliary << std::endl;

    double integral_sum_phi = 0.0;
    for (auto [phi_d_list, phi_d_dag_list] : phi_pair_list) {
      double integral_sum_n_left = 0.0;
      for (int n_left = 1; n_left < n; n_left++) {
        std::cout << " ------- tci start ------- " << std::endl;
        double integral_sum_iota   = 0.0;
        long count                 = 0;
        auto get_u_tau_max_element = [this, &count, &phi_d_list = phi_d_list, &phi_d_dag_list = phi_d_dag_list,
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

        // // find the pivot that gives a non-zero integrand
        // // this calculation can be performed only once and cached for all inchworm steps
        // std::cout << "searching for a non-zero integrand" << std::endl;
        // auto start_time_searching = std::chrono::high_resolution_clock::now();
        // std::vector<double> v_iota_s1{};
        // double u_tau_max_element_vs1 = 0;
        // int iota_pivot_index         = 0;
        // // parallelization is possible here
        // for (auto iota_pivot1 : all_iota_pivots) {
        //   std::vector<double> v_iota_s1_temp{};
        //   for (int i = 0; i < iota_pivot1.size(); i++) { v_iota_s1_temp.push_back(iotai[iota_pivot1[i]]); }
        //   for (int i = 0; i < v_pivot1.size(); i++) { v_iota_s1_temp.push_back(tp.vi[v_pivot1[i]]); }
        //   u_tau_max_element_vs1 = get_u_tau_max_element(v_iota_s1_temp);
        //   if (std::abs(u_tau_max_element_vs1) != 0) {
        //     v_iota_s1 = v_iota_s1_temp;
        //     break;
        //   }
        //   iota_pivot_index++;
        // }
        // auto end_time_searching            = std::chrono::high_resolution_clock::now();
        // auto duration_searching            = std::chrono::duration_cast<std::chrono::microseconds>(end_time_searching - start_time_searching).count();
        // auto duration_in_seconds_searching = static_cast<double>(duration_searching) / 1e6;
        // std::cout << "duration_in_seconds_searching: " << duration_in_seconds_searching << " seconds" << std::endl;
        // time_for_find_pivot += duration_in_seconds_searching;
        // std::cout << "finished searching for a non-zero integrand" << std::endl;
        // if (iota_pivot_index == all_iota_pivots.size()) {
        //   std::cout << "no non-zero integrand found, skip the integral" << std::endl;
        //   std::cout << " ------- tci finish ------- " << std::endl;
        //   continue;
        // }
        // int iota_pivot_index = 1;
        // auto pivot1           = all_iota_pivots[iota_pivot_index];
        // auto pivot1_to_append = v_pivot1;
        // pivot1.insert(pivot1.end(), pivot1_to_append.begin(), pivot1_to_append.end());
        auto pivot1_pre           = std::vector<int>(n, 1);
        auto pivot1_to_append_pre = v_pivot1_pre;
        pivot1_pre.insert(pivot1_pre.end(), pivot1_to_append_pre.begin(), pivot1_to_append_pre.end());
        // u_tau_max_element_vs1 = get_u_tau_max_element(v_iota_s1);
        // if (sp.debug > 1) {
        //   std::cout << std::endl;
        //   std::cout << "find a non-zero integrand" << std::endl;
        //   int mid1 = v_iota_s1.size() / 2;
        //   std::vector<double> iotas1(v_iota_s1.begin(), v_iota_s1.begin() + mid1);
        //   std::vector<double> vs1(v_iota_s1.begin() + mid1, v_iota_s1.end());
        //   std::vector<double> iota_d_list1     = get_elements(phi_d_list, iotas1);
        //   std::vector<double> iota_d_dag_list1 = get_elements(phi_d_dag_list, iotas1);
        //   std::vector<int> iota_d_list_int1(iota_d_list1.begin(), iota_d_list1.end());
        //   std::vector<int> iota_d_dag_list_int1(iota_d_dag_list1.begin(), iota_d_dag_list1.end());
        //   auto [taus_left1, taus_right1, taus1] = obtain_taus(vs1, n_left, sp.tau_split, sp.tau_max);
        //   print_pivot1(iota_d_list_int1, iota_d_dag_list_int1, get_elements(phi_d_list, taus1), get_elements(phi_d_dag_list, taus1),
        //                u_tau_max_element_vs1);
        // }

        long count_pre                 = 0;
        auto get_u_tau_max_element_pre = [this, &count_pre, &phi_d_list = phi_d_list, &phi_d_dag_list = phi_d_dag_list, &n_left,
                                          &auxi_height](const std::vector<double> &v_iota_s) {
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
          // double sin_term = linear_func_all( iotas, mp.n_phi);
          if (number_in_block_d != number_in_block_d_dag) { return auxi_height * sin_term; }
          //Note: if local Hamiltonian commute with density operators, further simplication could be applied: the operator of the same type can not be next to each other on the time axis
          auto [taus_left, taus_right, taus] = obtain_taus(vs, n_left, sp.tau_split, sp.tau_max);
          double integrand                   = evaluate_u_tau_max(sr.u_tau_max_zeroth_order, sp.tau_split, sp.tau_max, mp.all_d_ops, mp.all_d_dag_ops,
                                                                  mp.gf_block_shape, cp, mp.Delta_tau, mp.ad_imp, sr.u_interpolator, get_elements(phi_d_list, taus),
                                                                  get_elements(phi_d_dag_list, taus), iota_d_list, iota_d_dag_list, sp.bl_index, sp.subspace_index);
          count_pre++;
          double j = jacobian(taus_left, sp.tau_split, 0.0) * jacobian(taus_right, sp.tau_max, sp.tau_split);
          return auxi_height * sin_term + integrand * j;
        };

        //pretraining
        auto start_time_pre_training = std::chrono::high_resolution_clock::now();
        std::cout << "pretraining" << std::endl;
        std::cout << "iteration nEval LastSweepPivotError\n";
        auto ci_pre = xfac::CTensorCI2<double, double>(get_u_tau_max_element_pre, input_pre,
                                                       {.bondDim = tp.bond_dim, .reltol = reltol_test, .pivot1 = pivot1_auxiliary, .fullPiv = true});
        ci_pre.addPivots(ci_pre_auxiliary);
        double integral_pre = ci_pre.tt.sum(weight_pre);
        std::cout << "integral_pre (start): " << integral_pre << std::endl;
        std::cout << "integral_auxiliary: " << integral_auxiliary << std::endl;
        if (std::abs(integral_auxiliary - integral_pre) == 0) {
          std::cout << "pre_trained integral is too small, skip the integral" << std::endl;
          std::cout << "integral_auxiliary-integral_pre: " << integral_auxiliary - integral_pre << std::endl;
          std::cout << " ------- tci finish------- " << std::endl;
          continue;
        }
        // reset
        std::cout << "bond_dim: " << ci_pre.param.bondDim << std::endl;
        ci_pre                      = xfac::CTensorCI2<double, double>(get_u_tau_max_element_pre, input_pre,
                                                  {.bondDim = tp.bond_dim, .reltol = reltol_test, .pivot1 = pivot1_auxiliary, .fullPiv = false});
        int ci_count                = 0;
        double previous_pivot_error = -1E5;
        while (true) {
          ci_pre.iterate();
          // ci_pre.makeCanonical();
          auto last_pivot_error = ci_pre.pivotError[ci_pre.pivotError.size() - 1];
          std::cout << ci_count << " " << count_pre << " " << last_pivot_error << " " << std::endl;
          print_rank(ci_pre.tt);
          ci_count++;
          if (std::abs(last_pivot_error - previous_pivot_error) < tci_convergence_threshold && ci_count > 3) { break; }
          previous_pivot_error = last_pivot_error;
        }
        // ci_pre.makeCanonical();
        auto end_time_pre_training = std::chrono::high_resolution_clock::now();
        auto duration_pre_training = std::chrono::duration_cast<std::chrono::microseconds>(end_time_pre_training - start_time_pre_training).count();
        auto duration_in_seconds_pre_training = static_cast<double>(duration_pre_training) / 1e6;
        std::cout << "duration_in_seconds_pre_training: " << duration_in_seconds_pre_training << " seconds" << std::endl;
        time_for_pre_training += duration_in_seconds_pre_training;
        std::cout << "pretraining finished" << std::endl;

        integral_pre = ci_pre.tt.sum(weight_pre);
        std::cout << "integral_pre: " << integral_pre << std::endl;
        std::cout << "integral_auxiliary: " << integral_auxiliary << std::endl;
        double integral_diff = integral_pre - integral_auxiliary;
        std::cout << "integral_pre-integral_auxiliary: " << integral_diff << std::endl;
        if (integral_diff == 0) {
          std::cout << "pre_trained integral is too small, skip the integral" << std::endl;
          std::cout << "strange!!" << std::endl;
          std::cout << "integral_diff: " << integral_diff << std::endl;
          std::cout << " ------- tci finish------- " << std::endl;
          continue;
        }

        //training
        std::cout << "training" << std::endl;
        auto input_to_append = std::vector(n, tp.vi);
        auto input           = std::vector(n, iotai);
        input.insert(input.end(), input_to_append.begin(), input_to_append.end());
        auto weight_to_append = std::vector(n, tp.wi_v);
        auto weight           = std::vector(n, wi_iota);
        weight.insert(weight.end(), weight_to_append.begin(), weight_to_append.end());
        bool find_pivot1              = false;
        std::vector<int> pivot1_train = {};
        for (auto b = 0u; b < ci_pre.len() - 1; b++) {
          auto pivots = ci_pre.getPivotsAt(b);
          for (auto p : pivots) {
            auto mid = p.size() / 2;
            std::vector<double> v_iota_s(p.begin(), p.end() - mid);
            for (int i = 0; i < v_pivot1.size(); i++) { v_iota_s.push_back(tp.vi[int(tp.vi.size() / 2)]); }
            // print_vector(p);
            // print_vector(v_iota_s);
            auto val = get_u_tau_max_element(v_iota_s);
            if ( val != 0) {
              find_pivot1  = true;
              pivot1_train = p;
              std::cout << "found pivot1 with value: " << val << std::endl;
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
        std::cout << "iteration nEval LastSweepPivotError integral\n";
        auto ci = xfac::CTensorCI2<double, double>(get_u_tau_max_element, input,
                                                   {.bondDim = tp.bond_dim, .reltol = reltol_test, .pivot1 = pivot1_train, .fullPiv = false});
        print_rank(ci.tt);
        // ci_pre.iterate(2, 0);
        // ci_pre.iterate(2, 1);
        // ci.addPivots(ci_pre);
        for (auto b = 0u; b < ci.len() - 1; b++) {
          auto pivots = ci_pre.getPivotsAt(b);
          // auto first_half_pivots = std::vector(pivots.begin(), pivots.begin() + pivots.size() / 2);
          ci.myAddPivotsAt(pivots, b);
        }
        print_rank(ci.tt);
        // ci.makeCanonical();
        // print_rank(ci.tt);
        std::cout << "bond_dim: " << ci.param.bondDim << std::endl;
        double current_integral = 0.0;
        for (int i = 1; i <= tp.sweep_bound; i++) {
          ci.iterate();
          // ci.makeCanonical();
          current_integral      = ci.tt.sum(weight);
          auto last_pivot_error = ci.pivotError[ci.pivotError.size() - 1];
          std::cout << i << " " << count_pre << " " << last_pivot_error << " " << current_integral << std::endl;
          print_rank(ci.tt);
        }
        auto integral_element = current_integral;
        std::cout << " ------- tci finish ------- " << std::endl;

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
  std::cout << "all orders finished" << std::endl;
  std::cout << "time_for_find_pivot: " << time_for_find_pivot << " seconds" << std::endl;
  std::cout << "time_for_pre_training: " << time_for_pre_training << " seconds" << std::endl;
  std::cout << std::endl;
}
