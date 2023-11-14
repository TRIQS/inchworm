#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "./mode.hpp"

using namespace inchworm;

void ModeVertexFactorization::run_single_element() {

// TCI
  int n_phi                                 = std::accumulate(mp.gf_block_shape.begin(), mp.gf_block_shape.end(), 0);
  std::vector<double> iotai(n_phi);
  std::iota(iotai.begin(), iotai.end(), 0);
  auto wi_iota = std::vector(n_phi, 1.0);
  for (int order : sp.order_list) {
    auto start_time = std::chrono::high_resolution_clock::now();
    int n           = 2 * order;     // number of tau's
    std::vector<int> v_pivot1(n, 0); // pivots for tau, pivots for iota are added later
    std::vector<int> index_range(n);
    std::iota(index_range.begin(), index_range.end(), 0);
    auto phi_pair_list = get_all_phi(index_range); //gives all possible phi

    std::vector<int> iota_pivots(n, 0);
    std::vector<int> iota_pivots_range(n_phi);
    std::iota(iota_pivots_range.begin(), iota_pivots_range.end(), 0);
    std::vector<std::vector<int>> all_iota_pivots{};
    generate_combinations(iota_pivots_range, iota_pivots, 0, all_iota_pivots);
    // auto iota_pair_list     = get_all_iota(mp.gf_block_shape, order); //gives all possible iota
    double integral_sum_phi = 0.0;
    for (auto [phi_d_list, phi_d_dag_list] : phi_pair_list) {
      double integral_sum_n_left = 0.0;
      for (int n_left = 1; n_left < n; n_left++) {
        double integral_sum_iota = 0.0;
        // for (auto [iota_d_list, iota_d_dag_list] : iota_pair_list) {
        long count                 = 0;
        auto get_u_tau_max_element = [this,  &count, &phi_d_list = phi_d_list,
                                      &phi_d_dag_list = phi_d_dag_list, &n_left](const std::vector<double> &v_iota_s) {
          std::vector<double> vs{};
          std::vector<double> iotas{};
          vs.reserve(v_iota_s.size());
          iotas.reserve(v_iota_s.size());
          for (int i = 0; i < v_iota_s.size(); i++) {
            double intPart;
            double fracPart;
            fracPart = modf(v_iota_s[i], &intPart);
            vs.push_back(fracPart);
            iotas.push_back(intPart);
          }
          // std::cout << "vs: ";
          // print_vector(vs);
          // std::cout << "iotas: ";
          // print_vector(iotas);
          int mid_iota = iotas.size() / 2;
          // for(int i = 0 ; i < iotas.size(); i++){
          //   if(iotas[i]!=0){
          //     std::cout << "iota is not zero" << std::endl;
          //   }
          // }
          // int mid_iota = iotas.size() / 2;
          // std::vector<double> iota_d_list(iotas.begin(), iotas.begin() + mid_iota);
          // std::vector<double> iota_d_dag_list(iotas.begin() + mid_iota, iotas.end());
          std::vector<double> iota_d_list     = get_elements(phi_d_list, iotas);
          std::vector<double> iota_d_dag_list = get_elements(phi_d_dag_list, iotas);
          std::vector<int> iota_d_list_int(iota_d_list.begin(), iota_d_list.end());
          std::vector<int> iota_d_dag_list_int(iota_d_dag_list.begin(), iota_d_dag_list.end());
          std::vector<int> number_in_block_d     = generate_number_in_block(mp.gf_block_shape, iota_d_list_int);
          std::vector<int> number_in_block_d_dag = generate_number_in_block(mp.gf_block_shape, iota_d_dag_list_int);
          if (number_in_block_d != number_in_block_d_dag) { return 0.0; }

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
        double u_tau_max_element_vs1 = 0;

        // set pivot for iota
        int iota_pivot_index = 0;
        for (auto iota_pivot1 : all_iota_pivots) {
          std::vector<double> v_iota_s1_temp{};
          for (int i = 0; i < iota_pivot1.size(); i++) {
            double value = iota_pivot1[i] + tp.vi[v_pivot1[i]];
            v_iota_s1_temp.push_back(value);
          }
          u_tau_max_element_vs1 = get_u_tau_max_element(v_iota_s1_temp);
          if (u_tau_max_element_vs1 != 0) {
            v_iota_s1 = v_iota_s1_temp;
            break;
          }
          iota_pivot_index++;
        }
        if (iota_pivot_index == all_iota_pivots.size()) { continue; }

        std::vector<int> pivot1{};
        auto pivot1_to_append = all_iota_pivots[iota_pivot_index];
        pivot1.reserve(n);
        for (int i = 0; i < pivot1_to_append.size(); i++) { pivot1.push_back(v_pivot1[i] + pivot1_to_append[i] *tp.n_GK); }
        std::cout << "pivot1: ";
        print_vector(pivot1);
        std::cout << "v_iota_s1: ";
        print_vector(v_iota_s1);

        if (sp.debug) {
          std::vector<double> vs1{};
          std::vector<double> iotas1{};
          vs1.reserve(v_iota_s1.size());
          iotas1.reserve(v_iota_s1.size());
          for (int i = 0; i < v_iota_s1.size(); i++) {
            double intPart;
            double fracPart;
            fracPart = modf(v_iota_s1[i], &intPart);
            vs1.push_back(fracPart);
            iotas1.push_back(intPart);
          }
          int mid_iota1 = iotas1.size() / 2;
          std::vector<double> iota_d_list1(iotas1.begin(), iotas1.begin() + mid_iota1);
          std::vector<double> iota_d_dag_list1(iotas1.begin() + mid_iota1, iotas1.end());
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
          std::cout << "get_u_tau_max_element(pivot1): " << u_tau_max_element_vs1 << "\n" << std::endl;
        }
        if (u_tau_max_element_vs1 == 0) { continue; }

        double current_integral{0};
        double previous_integral{0};
        double integral_element{0};
        if (sp.debug) { std::cout << "iteration nEval LastSweepPivotError integral\n"; }
        std::vector<double> v_iota_i;
        std::vector<double> weight_i;
        for (int i = 0; i < n_phi; i++) {
          for (int j = 0; j < tp.vi.size(); j++) {
            v_iota_i.push_back(i + tp.vi[j]);
            weight_i.push_back( tp.wi_v[j]);
          }
        }

        std::vector<std::vector<double>> input  = std::vector(n, v_iota_i);
        std::vector<std::vector<double>> weight = std::vector(n, weight_i);
        // std::cout << "input: " << std::endl;
        // for (auto v : input) { print_vector(v); }
        // std::cout << "weight: " << std::endl;
        // for (auto v : weight) { print_vector(v); }
        if (tp.tci_prrlu) {
          auto ci = xfac::CTensorCI2<double, double>(get_u_tau_max_element, input, {.bond_dim = tp.bond_dim, .pivot1 = pivot1});
          for (int i = 0; i <  tp.sweep_bound; i++) {
            ci.iterate();
            // ci.makeCanonical();
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
            if (std::abs(current_integral - previous_integral) <  tp.error_bound && i > 5) { break; }
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
