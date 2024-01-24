#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "./mode.hpp"

using namespace inchworm;

void ModePartitionSin::run_single_element() {

  std::vector<double> iotai(mp.n_phi);
  std::iota(iotai.begin(), iotai.end(), 0);
  auto wi_iota = std::vector(mp.n_phi, 1.0);
  for (int order : sp.order_list) {
    auto max_weight                 = *std::max_element(tp.wi_v.begin(), tp.wi_v.end());
    double pre_integral_lower_bound = tp.integral_lower_bound / (std::pow(max_weight, order * 2) * std::pow(tp.wi_v.size(), order * 2));
    std::cout << "### order " << order << " ###" << std::endl;
    std::cout << "pre_integral_lower_bound: " << pre_integral_lower_bound << std::endl;
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
        double integral_sum_iota = 0.0;
        long count               = 0;

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

        std::vector<double> v_iota_s1{};
        double u_tau_max_element_vs1 = 0;

        // set pivot for iota
        int iota_pivot_index = 0;
        for (auto iota_pivot1 : all_iota_pivots) {
          std::vector<double> v_iota_s1_temp{};
          for (int i = 0; i < iota_pivot1.size(); i++) { v_iota_s1_temp.push_back(iotai[iota_pivot1[i]]); }
          for (int i = 0; i < v_pivot1.size(); i++) { v_iota_s1_temp.push_back(tp.vi[v_pivot1[i]]); }
          u_tau_max_element_vs1 = get_u_tau_max_element(v_iota_s1_temp);
          if (std::abs(u_tau_max_element_vs1) != 0) {
            v_iota_s1 = v_iota_s1_temp;
            break;
          }
          iota_pivot_index++;
        }
        if (iota_pivot_index == all_iota_pivots.size()) { continue; }

        auto pivot1           = all_iota_pivots[iota_pivot_index];
        auto pivot1_to_append = v_pivot1;
        pivot1.insert(pivot1.end(), pivot1_to_append.begin(), pivot1_to_append.end());

        u_tau_max_element_vs1 = get_u_tau_max_element(v_iota_s1);
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
        }
        if (u_tau_max_element_vs1 == 0) { continue; }

        double pre_factor = tp.auxi_height;
        double relto_test = tp.reltol;

        auto get_u_tau_max_element_pre = [this, &count, &phi_d_list = phi_d_list, &phi_d_dag_list = phi_d_dag_list, &n_left,
                                          &pre_factor](const std::vector<double> &v_iota_s) {
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
          double j        = jacobian(taus_left, sp.tau_split, 0.0) * jacobian(taus_right, sp.tau_max, sp.tau_split);
          double sin_term = sin_func_all({}, iotas, mp.n_phi);
          return pre_factor * sin_term + integrand * j;
        };

        auto get_u_tau_max_element_pre_abs = [this, &count, &phi_d_list = phi_d_list, &phi_d_dag_list = phi_d_dag_list, &n_left,
                                          &pre_factor](const std::vector<double> &v_iota_s) {
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
          double j        = jacobian(taus_left, sp.tau_split, 0.0) * jacobian(taus_right, sp.tau_max, sp.tau_split);
          double sin_term = sin_func_all({}, iotas, mp.n_phi);
          return std::abs(integrand * j);
        };

        //pretraining
        std::cout << "pretraining" << std::endl;
        auto input_to_append_pre  = std::vector(n, std::vector<double>{tp.vi[7]});
        auto input_pre            = std::vector(n, iotai);
        auto weight_artificial    = std::vector<double> {1.0};
        auto weight_to_append_pre = std::vector(n, weight_artificial);
        auto weight_pre           = std::vector(n, wi_iota);

        input_pre.insert(input_pre.end(), input_to_append_pre.begin(), input_to_append_pre.end());
        weight_pre.insert(weight_pre.end(), weight_to_append_pre.begin(), weight_to_append_pre.end());
        std::cout << "iteration nEval LastSweepPivotError\n";
        auto ci_pre = xfac::CTensorCI2<double, double>(get_u_tau_max_element_pre, input_pre,
                                                       {.bondDim = tp.bond_dim, .reltol = relto_test, .pivot1 = pivot1, .fullPiv = true});
        std::cout << "bond_dim: " << ci_pre.param.bondDim << std::endl;
        int ci_count                = 0;
        double previous_pivot_error = -1E5;
        int previous_pivot_count    = 0;
        while (true) {
          ci_pre.iterate();
          // ci_pre.makeCanonical();
          auto last_pivot_error = ci_pre.pivotError[ci_pre.pivotError.size() - 1];
          // auto last_pivot_error = ci_pre.trueError();
          std::cout << ci_count << " " << count << " " << last_pivot_error << " " << std::endl;
          print_rank(ci_pre.tt);
          // if (ci_count == 1 && last_pivot_error < tp.auxi_height) {
          //   std::cout << "probably too small, skip the integral" << std::endl;
          //   skip_integral = true;
          //   break;
          // }
          ci_count++;
          if (std::abs(last_pivot_error - previous_pivot_error) < 1e-20) { break; }
          if (ci_count == previous_pivot_count + 3) {
            previous_pivot_error = last_pivot_error;
            previous_pivot_count = ci_count;
          }
        }

        auto ci_pre_abs = xfac::CTensorCI2<double, double>(get_u_tau_max_element_pre_abs, input_pre,
                                                       {.bondDim = tp.bond_dim, .reltol = relto_test, .pivot1 = pivot1, .fullPiv = true});
        ci_pre_abs.addPivots(ci_pre);
        ci_pre_abs.makeCanonical();
        double pre_integral = ci_pre_abs.tt.sum(weight_pre);
        std::cout << "pre_integral: " << pre_integral << std::endl;
        if (std::abs(pre_integral) < pre_integral_lower_bound) {
          std::cout << "pre_trained integral is too small, skip the integral" << std::endl;
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

        std::cout << "iteration nEval LastSweepPivotError integral\n";
        auto ci = xfac::CTensorCI2<double, double>(get_u_tau_max_element, input,
                                                   {.bondDim = tp.bond_dim, .reltol = relto_test, .pivot1 = pivot1, .fullPiv = true});
        // for (auto b = 0u; b < ci.len() - 1; b++) { ci.myAddPivotsAt(ci_pre.getPivotsAt(b), b); }
        ci.addPivots(ci_pre);
        ci.makeCanonical();
        print_rank(ci.tt);
        std::cout << "bond_dim: " << ci.param.bondDim << std::endl;
        double current_integral = 0.0;
        for (int i = 1; i <= tp.sweep_bound; i++) {
          ci.iterate();
          ci.makeCanonical();
          current_integral      = ci.tt.sum(weight);
          auto last_pivot_error = ci.pivotError[ci.pivotError.size() - 1];
          std::cout << i << " " << count << " " << last_pivot_error << " " << current_integral << std::endl;
          print_rank(ci.tt);
        }
        auto integral_element = current_integral;

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
