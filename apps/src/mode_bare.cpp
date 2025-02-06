#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "mode.hpp"
#include "save.hpp"

using namespace inchworm;

void ModeBare::validate_input() {
  ModeBase::validate_input();
  if (gp.target != "propagator") {
    std::cerr << "bare mode: target must be propagator" << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

void ModeBare::print_summary() {
  if (sp.debug <= 0) return;
  ModeBase::print_summary();
  std::cout << "partition function exact: " << std::setw(10) << sr.partition_function_ref << std::endl;
  std::cout << "partition function exact * Z_imp_correction : " << std::setw(10) << sr.partition_function_ref * sr.Z_imp_correction << std::endl;
  std::cout << "partition function exact * Z_bath * Z_imp_correction * Z_bath_correction : " << std::setw(10)
            << sr.partition_function_ref * sr.Z_bath * sr.Z_imp_correction * sr.Z_bath_correction << std::endl;
  std::cout << std::left << std::setw(10) << "order" << std::setw(30) << "value" << std::setw(30) << "time(s)" << std::setw(30)
            << "time_find_pivot(s)" << std::setw(30) << "time_pretrain(s)" << std::setw(30) << "time_train(s)" << std::endl;
  auto partition_function_zeroth_order_ref = trace(sr.u_tau_zeroth_order);
  std::cout << std::setw(10) << "0" << std::setw(30) << partition_function_zeroth_order_ref << std::endl;
  for (int i = 0; i < sp.order_list.size(); i++) {
    std::cout << std::setw(10) << sp.order_list[i] << std::setw(30) << sr.integral_list[i][0] << std::setw(30) << (sr.statistics[0].time_order)[i]
              << std::endl;
  }
  double total_integral = 0.0;
  for (const auto &inner_vec : sr.integral_list) { total_integral += std::accumulate(inner_vec.begin(), inner_vec.end(), 0.0); }
  double sum_value = partition_function_zeroth_order_ref + total_integral;
  double sum_time = std::accumulate(sr.statistics[0].time_order.begin(), sr.statistics[0].time_order.end(), 0.0);
  std::cout << std::setw(10) << "sum:" << std::setw(30) << sum_value << std::setw(30) << sum_time << std::endl;
  std::cout << "sum * gp.Z_energy_shift_correction: " << std::setw(10) << sum_value * gp.Z_energy_shift_correction << std::endl;
  std::cout << "sum * Z_imp_correction * gp.Z_energy_shift_correction: " << std::setw(10)
            << sum_value * sr.Z_imp_correction * gp.Z_energy_shift_correction << std::endl;
  std::cout << "sum * Z_bath * Z_imp_correction * Z_bath_correction * gp.Z_energy_shift_correction: " << std::setw(10)
            << sum_value * sr.Z_bath * sr.Z_imp_correction * sr.Z_bath_correction * gp.Z_energy_shift_correction << std::endl;
}

void ModeBare::run() {
  std::cout << "### bare mode: start running ###" << std::endl;
  validate_input();
  ModeBase::prepare_eval();
  evaluate_propagator();
}

void ModeBare::evaluate_propagator() {
  sp.tau_split                 = 0.0;
  sp.tau_max                   = cp.beta;
  sp.bl_index                  = -1; // indicating doing trace instead of finding element
  sp.subspace_index            = -1;
  auto u_tau_zeroth_order_bare = make_bare_u_frame(mp.ad_imp, cp.beta, gp.exponent_u);
  sr.u_tau_zeroth_order        = u_tau_zeroth_order_bare;
  sp.use_bare_propagator       = true;
  // sr.u_interpolator    = interpolator_t<scalar_t>(sr.u_tau_ref, sr.u_tau_ref[0].mesh().size());
  ModeBase::evaluate();
}

void ModeBare::evaluate_greens_function() {}