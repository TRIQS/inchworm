#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "./mode.hpp"

using namespace inchworm;

void ModeTreeFactorization1::run_single_element() {

  std::vector<double> iotai(mp.n_phi);
  std::iota(iotai.begin(), iotai.end(), 0);
  auto wi_iota = std::vector(mp.n_phi, 1.0);
  for (int order : sp.order_list) {
    auto start_time = std::chrono::high_resolution_clock::now();
    int n           = 2 * order; // number of tau's
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
        auto get_u_tau_max_element = [this, &count, &phi_d_list = phi_d_list, &phi_d_dag_list = phi_d_dag_list, &n_left,
                                      &n](const std::vector<double> &v_iota_s) {
          std::vector<double> vs{};
          std::vector<double> iotas{};
          for (int ind_tau = 0; ind_tau < n; ++ind_tau) {
            bool found       = false;
            bool found_twice = false;
            for (int ind_phi = 0; ind_phi < mp.n_phi; ++ind_phi) {
              if (v_iota_s[ind_tau * mp.n_phi + ind_phi] > 0) {
                vs.push_back(v_iota_s[ind_tau * mp.n_phi + ind_phi]);
                iotas.push_back(ind_phi);
                if (found) { found_twice = true; }
                found = true;
              }
            }
            if (!found || found_twice) { return 0.0; }
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
        double u_tau_max_element_vs1 = 0;
        int iota_pivot_index         = 0;
        bool found_pivot1            = false;
        std::vector<int> pivot1_init(mp.n_phi * n, tp.n_GK); // set the pivot to be indicator of not selected
        std::vector<double> v_iota_s1_init(mp.n_phi * n, -1.0);
        std::vector<int> pivot1 = pivot1_init;
        std::vector<double> v_iota_s1 = v_iota_s1_init;

        for (auto iota_pivot1 : all_iota_pivots) {
          auto pivot1_temp    = pivot1_init;
          auto v_iota_s1_temp = v_iota_s1_init;
          for (int i = 0; i < iota_pivot1.size(); i++) {
            pivot1_temp[i * mp.n_phi + iota_pivot1[i]]    = 0;
            v_iota_s1_temp[i * mp.n_phi + iota_pivot1[i]] = tp.vi[0];
          }
          double u_tau_max_element_vs1_temp = get_u_tau_max_element(v_iota_s1_temp);
          if (u_tau_max_element_vs1_temp > u_tau_max_element_vs1) {
            u_tau_max_element_vs1 = u_tau_max_element_vs1_temp;
            pivot1                = pivot1_temp;
            v_iota_s1             = v_iota_s1_temp;
            found_pivot1          = true;
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

        if (sp.debug > 1) {
          std::vector<double> vs1{};
          std::vector<double> iotas1{};
          for (int ind_tau = 0; ind_tau < tp.n_GK; ++ind_tau) {
            for (int ind_phi = 0; ind_phi < mp.n_phi; ++ind_phi) {
              if (v_iota_s1[ind_tau * mp.n_phi + ind_phi] > 0) {
                vs1.push_back(v_iota_s1[ind_tau * mp.n_phi + ind_phi]);
                iotas1.push_back(ind_phi);
                break;
              }
            }
          }
          std::cout << std::endl;
          std::vector<double> iota_d_list1     = get_elements(phi_d_list, iotas1);
          std::vector<double> iota_d_dag_list1 = get_elements(phi_d_dag_list, iotas1);
          std::vector<int> iota_d_list_int1(iota_d_list1.begin(), iota_d_list1.end());
          std::vector<int> iota_d_dag_list_int1(iota_d_dag_list1.begin(), iota_d_dag_list1.end());
          auto [taus_left1, taus_right1, taus1] = obtain_taus(vs1, n_left, sp.tau_split, sp.tau_max);
          print_pivot1(iota_d_list_int1, iota_d_dag_list_int1, get_elements(phi_d_list, taus1), get_elements(phi_d_dag_list, taus1),
                       u_tau_max_element_vs1);
        }
        if (u_tau_max_element_vs1 == 0) { continue; }

        std::vector<double> v_iota_i = tp.vi;
        std::vector<double> weight_i = tp.wi_v;
        v_iota_i.push_back(-1.0);
        weight_i.push_back(1.0);

        //for debuging
        if (order == 2) {
          std::cout << "------------" << std::endl;
          std::cout << "n_left " << n_left << std::endl;
          std::cout << "phi_d_list: " << std::endl;
          for (auto i : phi_d_list) { std::cout << i << " "; }
          std::cout << std::endl;
          std::cout << "phi_d_dag_list: " << std::endl;
          for (auto i : phi_d_dag_list) { std::cout << i << " "; }
          std::cout << std::endl;
          std::vector<double> v_iota_s2{0 + tp.vi[0], 0 + tp.vi[0], 0 + tp.vi[0], 0 + tp.vi[0]};
          std::cout << "get_u_tau_max_element(0,0,0,0): " << get_u_tau_max_element(v_iota_s2) << std::endl;
          std::vector<double> v_iota_s3{0 + tp.vi[0], 0 + tp.vi[0], 0 + tp.vi[0], 1 + tp.vi[0]};
          std::cout << "get_u_tau_max_element(0,0,0,1): " << get_u_tau_max_element(v_iota_s3) << std::endl;
          std::vector<double> v_iota_s4{0 + tp.vi[0], 0 + tp.vi[0], 1 + tp.vi[0], 0 + tp.vi[0]};
          std::cout << "get_u_tau_max_element(0,0,1,0): " << get_u_tau_max_element(v_iota_s4) << std::endl;
          std::vector<double> v_iota_s5{0 + tp.vi[0], 0 + tp.vi[0], 1 + tp.vi[0], 1 + tp.vi[0]};
          std::cout << "get_u_tau_max_element(0,0,1,1): " << get_u_tau_max_element(v_iota_s5) << std::endl;
          std::vector<double> v_iota_s6{1 + tp.vi[0], 1 + tp.vi[0], 0 + tp.vi[0], 0 + tp.vi[0]};
          std::cout << "get_u_tau_max_element(1,1,0,0): " << get_u_tau_max_element(v_iota_s6) << std::endl;
          std::vector<double> v_iota_s7{0 + tp.vi[0], 1 + tp.vi[0], 0 + tp.vi[0], 1 + tp.vi[0]};
          std::cout << "get_u_tau_max_element(0,1,0,1): " << get_u_tau_max_element(v_iota_s7) << std::endl;
          std::cout << "------------" << std::endl;
        }

        std::vector<std::vector<double>> input  = std::vector(n * mp.n_phi, v_iota_i);
        std::vector<std::vector<double>> weight = std::vector(n * mp.n_phi, weight_i);
        double integral_element                 = do_TCI<double, double>(get_u_tau_max_element, input, weight, pivot1, tp.sweep_bound, tp.bond_dim,
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
