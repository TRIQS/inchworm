#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "./mode.hpp"

using namespace inchworm;

void ModeBare::validate_input() {
  ModeBase::validate_input();
  if (gp.target != "propagator") {
    std::cerr << "bare mode: target must be propagator" << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

void ModeBare::run() {
  std::cout << "### bare mode: start running ###" << std::endl;
  validate_input();
  evaluate_propagator();
}

void ModeBare::evaluate_propagator() {
    sp.tau_split      = 0.0;
    sp.tau_max        = cp.beta;
    sp.bl_index       = -1; // indicating doing trace instead of finding element
    sp.subspace_index = -1;
    sr.u_tau_zeroth_order = sr.u_tau_zeroth_order_bare;
    sp.use_bare_propagator = true;
    // sr.u_interpolator    = interpolator_t<scalar_t>(sr.u_tau_ref, sr.u_tau_ref[0].mesh().size());
    ModeBase::evaluate();
}

// void ModeBare::evaluate_propagator() {
//   //In the bare mode, the propagator refers to partition function
//   //discrete index that can be summed or be in tci: phi, iota (at most 2 loops for summation)
//   auto loop1 = Loop("empty", std::vector<int>{0});
//   auto loop2 = Loop("empty", std::vector<int>{0});
//   cv_func change_variable;
//   jb_func jacobian;
//   sp.tau_split      = 0.0;
//   sp.tau_max        = cp.beta;
//   sp.bl_index       = -1; // indicating doing trace instead of finding element
//   sp.subspace_index = -1;
//   int n_left        = 0;
//   if (tp.mapping_v == 0) {
//     change_variable = change_variable0;
//     jacobian        = jacobian0;
//   } else if (tp.mapping_v == 1) {
//     change_variable = change_variable1;
//     jacobian        = jacobian1;
//   }
//   else if (tp.mapping_v == 2) {
//     change_variable = change_variable2;
//     jacobian        = jacobian2;
//   }
//   else if (tp.mapping_v == 3) {
//     change_variable = change_variable3;
//     jacobian        = jacobian3;
//   }
//   else {
//     std::cerr << "invalid mapping_v" << std::endl;
//     std::exit(EXIT_FAILURE);
//   }
//   for (int order : sp.order_list) {
//     if (sp.debug) std::cout << "order: " << order << std::endl;
//     int n = 2 * order; // number of operators
//     // generate valid phi and iota pairs
//     std::vector<int> index_range(n);
//     std::iota(index_range.begin(), index_range.end(), 0);
//     auto phi_pair_list = get_all_phi(index_range);
//     std::vector<int> phi_list(phi_pair_list.size());
//     std::iota(phi_list.begin(), phi_list.end(), 0);
//     auto iota_pair_list = get_all_iota(mp.gf_block_shape, order);
//     std::vector<int> iota_list(iota_pair_list.size());
//     std::iota(iota_list.begin(), iota_list.end(), 0);
//     // generate valid iota pivots (used when the integral_variable == v_iota)
//     std::vector<int> iota_pivots(n, 0); // this is an intermediate variable for generating all possible iota
//     std::vector<int> iota_pivots_range(mp.n_phi);
//     std::iota(iota_pivots_range.begin(), iota_pivots_range.end(), 0);
//     std::vector<std::vector<int>> all_iota_pivots{};
//     generate_combinations(iota_pivots_range, iota_pivots, 0, all_iota_pivots);
//     if (gp.integrand == "plain" && gp.integral_variable == "v" && gp.tci_shape == "plain") {
//       loop1 = Loop("phi", phi_list);
//       loop2 = Loop("iota", iota_list);
//     } else if (gp.integrand == "sum_phi" && gp.integral_variable == "v_iota" && gp.tci_shape == "vertex") {
//     } else {
//       std::cerr << "bare mode: invalid integrand or integral_variable or tci_shape" << std::endl;
//       std::exit(EXIT_FAILURE);
//     }
//     auto start_time        = std::chrono::high_resolution_clock::now();
//     double time_find_pivot = 0.0;
//     double time_pretrain   = 0.0;
//     double time_train      = 0.0;
//     loop1.value            = 0;

//     for (auto val1 : loop1.container) {
//       loop2.value = 0;
//       for (auto val2 : loop2.container) {
//         std::vector<int> phi_d_list{};
//         std::vector<int> phi_d_dag_list{};
//         std::vector<int> iota_d_list{};
//         std::vector<int> iota_d_dag_list{};
//         int id_phi  = -1;
//         int id_iota = -1;
//         if (loop1.name == "phi")
//           id_phi = val1;
//         else if (loop2.name == "phi")
//           id_phi = val2;
//         if (loop1.name == "iota")
//           id_iota = val1;
//         else if (loop2.name == "iota")
//           id_iota = val2;
//         if (id_phi != -1) {
//           std::tie(phi_d_list, phi_d_dag_list) = phi_pair_list[id_phi];
//         } else if (gp.integrand == "sum_phi") {
//         } else {
//           std::cerr << "invalid combination of integrand, integral_variable, and tci_shape" << std::endl;
//           std::exit(EXIT_FAILURE);
//         }
//         if (id_iota != -1) {
//           std::tie(iota_d_list, iota_d_dag_list) = iota_pair_list[id_iota];
//         } else if (gp.integral_variable == "v_iota") {
//         } else {
//           std::cerr << "invalid combination of integrand, integral_variable, and tci_shape" << std::endl;
//           std::exit(EXIT_FAILURE);
//         }

//         long count     = 0;
//         auto integrand = [this, &count, &phi_d_list = phi_d_list, &phi_d_dag_list = phi_d_dag_list, &iota_d_list = iota_d_list,
//                           &iota_d_dag_list = iota_d_dag_list, &n_left, &change_variable, &jacobian,
//                           &phi_pair_list](const std::vector<double> &variables) -> double {
//           std::vector<double> taus_left{};
//           std::vector<double> taus_right{};
//           std::vector<double> taus{};
//           count++;
//           double integrand = 0.0;
//           if (gp.integral_variable == "v" && gp.integrand == "plain") {
//             std::tie(taus_left, taus_right, taus) = obtain_taus(variables, n_left, sp.tau_split, sp.tau_max, change_variable);
//             integrand = evaluate_u_tau_max(sr.u_tau_zeroth_order_bare, sp.tau_split, sp.tau_max, mp.all_d_ops, mp.all_d_dag_ops, mp.gf_block_shape,
//                                            cp, mp.Delta_tau, mp.ad_imp, sr.u_interpolator_ref, get_elements(phi_d_list, taus),
//                                            get_elements(phi_d_dag_list, taus), iota_d_list, iota_d_dag_list, sp.bl_index, sp.subspace_index);
//             double j  = jacobian(taus_right, sp.tau_max, sp.tau_split);
//             return integrand * j;
//           } else if (gp.integral_variable == "v_iota" && gp.integrand == "sum_phi" && gp.tci_shape == "vertex") {
//             std::vector<double> vs{};
//             std::vector<double> iotas{};
//             vs.reserve(variables.size());
//             iotas.reserve(variables.size());
//             for (int i = 0; i < variables.size(); i++) {
//               double int_part;
//               double frac_part;
//               frac_part = modf(variables[i], &int_part);
//               vs.push_back(frac_part);
//               iotas.push_back(int_part);
//             }
//             std::tie(taus_left, taus_right, taus) = obtain_taus(vs, n_left, sp.tau_split, sp.tau_max, change_variable);
//             for (auto [phi_d_list, phi_d_dag_list] : phi_pair_list) {
//               std::vector<double> iota_d_list     = get_elements(phi_d_list, iotas);
//               std::vector<double> iota_d_dag_list = get_elements(phi_d_dag_list, iotas);
//               std::vector<int> iota_d_list_int(iota_d_list.begin(), iota_d_list.end());
//               std::vector<int> iota_d_dag_list_int(iota_d_dag_list.begin(), iota_d_dag_list.end());
//               std::vector<int> number_in_block_d     = generate_number_in_block(mp.gf_block_shape, iota_d_list_int);
//               std::vector<int> number_in_block_d_dag = generate_number_in_block(mp.gf_block_shape, iota_d_dag_list_int);
//               if (number_in_block_d != number_in_block_d_dag) { continue; }
//               auto integrand_phi =
//                  evaluate_u_tau_max(sr.u_tau_zeroth_order_bare, sp.tau_split, sp.tau_max, mp.all_d_ops, mp.all_d_dag_ops, mp.gf_block_shape, cp,
//                                     mp.Delta_tau, mp.ad_imp, sr.u_interpolator_ref, get_elements(phi_d_list, taus),
//                                     get_elements(phi_d_dag_list, taus), iota_d_list, iota_d_dag_list, sp.bl_index, sp.subspace_index);
//               double j = jacobian(taus_right, sp.tau_max, sp.tau_split);
//               integrand += integrand_phi * j;
//             }
//             return integrand;
//           } else {
//             std::cerr << "not implemented" << std::endl;
//             std::exit(EXIT_FAILURE);
//           }
//         };
//         //training
//         std::vector<std::vector<double>> input{};
//         std::vector<std::vector<double>> weight{};
//         std::vector<int> init_pivot{};
//         if (gp.integral_variable == "v" && gp.integrand == "plain" && gp.tci_shape == "plain") {
//           for (int i = 0; i < n; i++) {
//             input.push_back(tp.v_value);
//             weight.push_back(tp.v_weight);
//             init_pivot.push_back(0);
//           }
//         } else if (gp.integral_variable == "v_iota" && gp.integrand == "sum_phi" && gp.tci_shape == "vertex") {
//           std::vector<double> v_iota_value;
//           std::vector<double> v_iota_weight;
//           for (int i = 0; i < mp.n_phi; i++) {
//             for (int j = 0; j < tp.v_value.size(); j++) {
//               v_iota_value.push_back(i + tp.v_value[j]);
//               v_iota_weight.push_back(tp.v_weight[j]);
//             }
//           }
//           for (int i = 0; i < n; i++) {
//             input.push_back(v_iota_value);
//             weight.push_back(v_iota_weight);
//             init_pivot.push_back(0);
//           }
//         } else {
//           std::cerr << "not implemented" << std::endl;
//           std::exit(EXIT_FAILURE);
//         }

//         std::vector<double> init_input{};
//         for (int i = 0; i < init_pivot.size(); i++) { init_input.push_back(input[i][init_pivot[i]]); }
//         double init_integrand = integrand(init_input);
//         if (init_integrand == 0) { continue; }
//         double integral = do_TCI<double, double>(integrand, input, weight, init_pivot, count, tp.sweep_bound, tp.bond_dim, tp.reltol, tp.fullPiv,
//                                                  tp.tci_prrlu, tp.error_type, tp.error_eval, tp.convergence_bound, tp.convergence_iter, sp.debug);
//         loop2.value += integral;
//       } // end of loop2
//       loop1.value += loop2.value;
//     } // end of loop1
//     sr.integral_list.push_back(loop1.value);
//     auto end_time = std::chrono::high_resolution_clock::now();
//     sr.calculation_time_list.push_back(std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() / 1e6);
//     sr.find_pivot_time_list.push_back(time_find_pivot);
//     sr.pretrain_time_list.push_back(time_pretrain);
//     sr.train_time_list.push_back(time_train);
//   } // end of order loop
//   std::cout << "completed" << std::endl;
// }

void ModeBare::evaluate_greens_function() {}