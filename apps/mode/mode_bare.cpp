#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "./mode.hpp"

using namespace inchworm;

void ModeBare::validate_input() {
  ModeBase::validate_input();
  if (gp.target != "propagator" || gp.model_type != 0) {
    std::cerr << "bare mode: target must be propagator and model_type must be 0 (discerete bath)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

void ModeBare::run() {
  std::cout << "### bare mode: start running ###" << std::endl;
  validate_input();
  evaluate_propagator();
}

void ModeBare::evaluate_propagator() {
  //In the bare mode, the propagator refers to partition function
  //discrete index that can be summed or be in tci: phi, iota (at most 2 loops for summation)
  auto loop1 = Loop("empty", std::vector<int>{0});
  auto loop2 = Loop("empty", std::vector<int>{0});
  cv_func change_variable;
  jb_func jacobian;
  sp.tau_split      = 0.0;
  sp.tau_max        = cp.beta;
  sp.bl_index       = -1; // indicating doing trace instead of finding element
  sp.subspace_index = -1;
  int n_left        = 0;
  if (tp.mapping_v == 0) {
    change_variable = change_variable0;
    jacobian        = jacobian0;
  } else if (tp.mapping_v == 1) {
    change_variable = change_variable1;
    jacobian        = jacobian1;
  } else {
    std::cerr << "invalid mapping_v" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  for (int order : sp.order_list) {
    int n = 2 * order; // number of operators
    std::vector<int> index_range(n);
    std::iota(index_range.begin(), index_range.end(), 0);
    auto phi_pair_list = get_all_phi(index_range);
    std::vector<int> phi_list(phi_pair_list.size());
    std::iota(phi_list.begin(), phi_list.end(), 0);
    auto iota_pair_list = get_all_iota(mp.gf_block_shape, order);
    std::vector<int> iota_list(iota_pair_list.size());
    std::iota(iota_list.begin(), iota_list.end(), 0);
    if (gp.integrand == "plain" && gp.integral_variable == "v" && gp.tci_shape == "plain") {
      loop1 = Loop("phi", phi_list);
      loop2 = Loop("iota", iota_list);
    } else if (gp.integrand == "sum_phi" && gp.integral_variable == "v_iota" && gp.tci_shape == "vertex") {
    } else {
      std::cerr << "bare mode: invalid integrand or integral_variable or tci_shape" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    auto start_time        = std::chrono::high_resolution_clock::now();
    double time_find_pivot = 0.0;
    double time_pretrain   = 0.0;
    double time_train      = 0.0;
    loop1.value            = 0;

    for (auto val1 : loop1.container) {
      loop2.value = 0;
      for (auto val2 : loop2.container) {
        std::vector<int> phi_d_list{};
        std::vector<int> phi_d_dag_list{};
        std::vector<int> iota_d_list{};
        std::vector<int> iota_d_dag_list{};
        int id_phi  = -1;
        int id_iota = -1;
        if (loop1.name == "phi")
          id_phi = val1;
        else if (loop2.name == "phi")
          id_phi = val2;
        if (loop1.name == "iota")
          id_iota = val1;
        else if (loop2.name == "iota")
          id_iota = val2;
        if (id_phi != -1) {
          std::tie(phi_d_list, phi_d_dag_list) = phi_pair_list[id_phi];
        } else if (gp.integrand == "sum_phi") {
        } else {
          std::cerr << "invalid combination of integrand, integral_variable, and tci_shape" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (id_iota != -1) {
          std::tie(iota_d_list, iota_d_dag_list) = iota_pair_list[id_iota];
        } else if (gp.integral_variable == "v_iota") {
        } else {
          std::cerr << "invalid combination of integrand, integral_variable, and tci_shape" << std::endl;
          std::exit(EXIT_FAILURE);
        }

        long count     = 0;
        auto integrand = [this, &count, &phi_d_list = phi_d_list, &phi_d_dag_list = phi_d_dag_list, &iota_d_list = iota_d_list,
                          &iota_d_dag_list = iota_d_dag_list, &n_left, &change_variable, &jacobian](const std::vector<double> &variables) -> double {
          std::vector<double> taus_left{};
          std::vector<double> taus_right{};
          std::vector<double> taus{};
          double integrand = 0.0;
          if (gp.integral_variable == "v" && gp.integrand == "plain") {
            std::tie(taus_left, taus_right, taus) = obtain_taus(variables, n_left, sp.tau_split, sp.tau_max, change_variable);
            integrand = evaluate_u_tau_max(sr.u_tau_zeroth_order_ref, sp.tau_split, sp.tau_max, mp.all_d_ops, mp.all_d_dag_ops, mp.gf_block_shape, cp,
                                           mp.Delta_tau, mp.ad_imp, sr.u_interpolator_ref, get_elements(phi_d_list, taus),
                                           get_elements(phi_d_dag_list, taus), iota_d_list, iota_d_dag_list, sp.bl_index, sp.subspace_index);
          } else if (gp.integral_variable == "v_iota" && gp.integrand == "sum_phi" && gp.tci_shape == "vertex") {
            std::cerr << "not implemented" << std::endl;
            std::exit(EXIT_FAILURE);
          } else {
            std::cerr << "not implemented" << std::endl;
            std::exit(EXIT_FAILURE);
          }

          count++;
          double j = jacobian(taus_right, sp.tau_max, sp.tau_split);
          return integrand * j;
        };
        //training
        std::vector<std::vector<double>> input{};
        std::vector<std::vector<double>> weight{};
        std::vector<int> init_pivot{};
        if (gp.integral_variable == "v") {
          for (int i = 0; i < n; i++) {
            input.push_back(tp.v_value);
            weight.push_back(tp.v_weight);
            init_pivot.push_back(0);
          }
        } else {
          std::cerr << "not implemented" << std::endl;
          std::exit(EXIT_FAILURE);
        }

        std::vector<double> init_input{};
        for (int i = 0; i < init_pivot.size(); i++) { init_input.push_back(input[i][init_pivot[i]]); }
        double init_integrand = integrand(init_input);
        if (init_integrand == 0) { continue; }
        double integral = do_TCI<double, double>(integrand, input, weight, init_pivot, count, tp.sweep_bound, tp.bond_dim, tp.reltol, tp.fullPiv,
                                                 tp.tci_prrlu, tp.error_type, tp.error_eval, tp.convergence_bound, tp.convergence_iter, sp.debug);
        loop2.value += integral;
      } // end of loop2
      loop1.value += loop2.value;
    } // end of loop1
    sr.integral_list.push_back(loop1.value);
    auto end_time = std::chrono::high_resolution_clock::now();
    sr.calculation_time_list.push_back(std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() / 1e6);
    sr.find_pivot_time_list.push_back(time_find_pivot);
    sr.pretrain_time_list.push_back(time_pretrain);
    sr.train_time_list.push_back(time_train);
  } // end of order loop
    std::cout << "completed" << std::endl;
}

void ModeBare::evaluate_greens_function() {}