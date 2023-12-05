#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "./mode.hpp"

using namespace inchworm;

void ModeVertexFactorizationPivots::run_single_element() {

  std::vector<double> iotai(mp.n_phi);
  std::iota(iotai.begin(), iotai.end(), 0);
  auto wi_iota = std::vector(mp.n_phi, 1.0);
  for (int order : sp.order_list) {
    auto start_time = std::chrono::high_resolution_clock::now();
    int n           = 2 * order;     // number of tau's
    std::vector<int> v_pivot1(n, 7); // pivots for tau, pivots for iota are added later
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
          vs.reserve(v_iota_s.size());
          iotas.reserve(v_iota_s.size());
          for (int i = 0; i < v_iota_s.size(); i++) {
            double int_part;
            double frac_part;
            frac_part = modf(v_iota_s[i], &int_part);
            vs.push_back(frac_part);
            iotas.push_back(int_part);
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

          // double sum_iota_d = std::accumulate(iota_d_list.begin(), iota_d_list.end(), 0.0);
          // double sum_iota_d_dag = std::accumulate(iota_d_dag_list.begin(), iota_d_dag_list.end(), 0.0);
          // if(sum_iota_d != sum_iota_d_dag && integrand !=0){
          //   std::cout << "strange" << std::endl;
          // }

          double j = jacobian(taus_left, sp.tau_split, 0.0) * jacobian(taus_right, sp.tau_max, sp.tau_split);
          return integrand * j;
        };

        std::vector<double> v_iota_s1;
        double u_tau_max_element_vs1 = 0;

        // set pivot for iota
        int iota_pivot_index = 0;
        bool found_pivot1    = false;
        std::vector<int> valid_iota_index{};
        std::vector<double> valid_iota_value{};
        for (auto iota_pivot1 : all_iota_pivots) {
          std::vector<double> v_iota_s1_temp{};
          for (int i = 0; i < iota_pivot1.size(); i++) {
            double value = iota_pivot1[i] + tp.vi[v_pivot1[i]];
            v_iota_s1_temp.push_back(value);
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
        // for debugging, print the size of valid_iota_index
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
          std::vector<int> pivot1_temp{};
          auto pivot1_to_append = all_iota_pivots[i];
          pivot1_temp.reserve(n);
          for (int j = 0; j < pivot1_to_append.size(); j++) { pivot1_temp.push_back(v_pivot1[j] + pivot1_to_append[j] * tp.n_GK); }
          valid_pivots.push_back(pivot1_temp);
        }
        auto pivot1 = valid_pivots[0];

        if (sp.debug > 1) {
          std::vector<double> vs1{};
          std::vector<double> iotas1{};
          vs1.reserve(v_iota_s1.size());
          iotas1.reserve(v_iota_s1.size());
          for (int i = 0; i < v_iota_s1.size(); i++) {
            double int_part;
            double frac_part;
            frac_part = modf(v_iota_s1[i], &int_part);
            vs1.push_back(frac_part);
            iotas1.push_back(int_part);
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

        std::vector<double> v_iota_i;
        std::vector<double> weight_i;
        for (int i = 0; i < mp.n_phi; i++) {
          for (int j = 0; j < tp.vi.size(); j++) {
            v_iota_i.push_back(i + tp.vi[j]);
            weight_i.push_back(tp.wi_v[j]);
          }
        }

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

        std::vector<std::vector<double>> input  = std::vector(n, v_iota_i);
        std::vector<std::vector<double>> weight = std::vector(n, weight_i);
        double integral_element =
           do_TCI_add_pivots<double, double>(get_u_tau_max_element, input, weight, pivot1, tp.sweep_bound, tp.bond_dim, tp.integral_error_bound,
                                             tp.pivot_error_bound, tp.tci_prrlu, sp.debug, count, valid_pivots);
        integral_sum_iota += integral_element;
        integral_sum_n_left += integral_sum_iota;
        //for debugging
        // for debugging
        // if (n_left == 2 && phi_d_list[0] == 1 && phi_d_list[1] == 3 && phi_d_dag_list[0] == 0 && phi_d_dag_list[1] == 2) {
        //   auto ci = xfac::CTensorCI2<double, double>(get_u_tau_max_element, input,
        //                                              {.bond_dim = tp.bond_dim, .reltol = 1e-18, .do_full_search = true, .pivot1 = pivot1});
        //   std::cout << "bond_dim: " << ci.param.bond_dim << std::endl;
        //   double current_integral = 0;
        //   double last_pivot_error = 0;
        //   for (int i = 1; i <= tp.sweep_bound; i++) {
        //     ci.iterate();
        //     current_integral = ci.tt.sum(weight);
        //     last_pivot_error = ci.pivotError[ci.pivotError.size() - 1];
        //     std::cout << i << " " << count << " " << last_pivot_error << " " << current_integral << std::endl;
        //   }
        //   std::ofstream outfile("./tci_plain.txt");
        //   std::cout << "writing tci" << std::endl;
        //   for (int id0 = 0; id0 < v_iota_i.size(); id0++) {
        //     for (int id1 = 0; id1 < v_iota_i.size(); id1++) {
        //       for (int id2 = 0; id2 < v_iota_i.size(); id2++) {
        //         for (int id3 = 0; id3 < v_iota_i.size(); id3++) {
        //           // double element = get_u_tau_max_element({v_iota_i[id0], v_iota_i[id1], v_iota_i[id2], v_iota_i[id3]});
        //           double element = ci.tt.eval({id0, id1, id2, id3});
        //           outfile << element << " ";
        //         }
        //       }
        //     }
        //   }
        //   outfile.close();
        // }
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
