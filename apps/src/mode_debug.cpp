#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "mode.hpp"
#include "save.hpp"

using namespace inchworm;

void ModeDebug::validate_input() {

  ModeBase::validate_input();
  if (gp.model_type != 0) {
    std::cerr << "debug mode: model_type must be 0 (discerete bath)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (sp.tau_max < 0 || sp.tau_max > cp.beta || sp.tau_split < 0 || sp.tau_split > sp.tau_max) {
    std::cerr << "debug mode: invalid tau_max or tau_split" << std::endl;
    std::cerr << "tau_max: " << sp.tau_max << std::endl;
    std::cout << "beta: " << cp.beta << std::endl;
    std::cerr << "tau_split: " << sp.tau_split << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (sp.bl_index < 0 || sp.bl_index >= sr.u_tau_ref.size()) {
    std::cerr << "debug mode: invalid bl_index" << std::endl;
    std::cerr << "sp.bl_index: " << sp.bl_index << std::endl;
    std::cerr << "It should be in the range of [0, " << sr.u_tau_ref.size() << ")" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (sp.subspace_index < 0 || sp.subspace_index >= sr.u_tau_ref[sp.bl_index].target_shape()[0] * sr.u_tau_ref[sp.bl_index].target_shape()[1]) {
    std::cerr << "debug mode: invalid subspace_index" << std::endl;
    std::cerr << "sp.subspace_index: " << sp.subspace_index << std::endl;
    std::cerr << "It should be in the range of [0, " << sr.u_tau_ref[sp.bl_index].target_shape()[0] * sr.u_tau_ref[sp.bl_index].target_shape()[1]
              << ")" << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

void ModeDebug::print_summary() {
  if (sp.debug <= 0) return;
  ModeBase::print_params();
  if (gp.target == "propagator") {

    auto [bl, i, j] = bl2_to_bl3(sp.bl_index, sp.subspace_index, mp.ad_imp.get_subspace_dims());
    std::cout << "#### debug results ####" << std::endl;
    std::cout << "---- order contribution ----" << std::endl;
    std::cout << std::left << std::setw(10) << "order" << std::setw(30) << "value" << std::setw(30) << "time(s)" << std::setw(30) << std::endl;
    std::cout << std::setw(10) << "0" << std::setw(30) << sr.u_tau_zeroth_order[sp.bl_index](i, j) << std::endl;
    for (int i = 0; i < sp.order_list.size(); i++) {
      std::cout << std::setw(10) << sp.order_list[i] << std::setw(30) << sr.integral_list[i][0] << std::setw(30) << (sr.statistics[0].time_order)[i]
                << std::endl;
    }
    double total_integral        = 0.0;
    gp.Z_energy_shift_correction = std::exp(gp.exponent_u * cp.beta);
    for (const auto &inner_vec : sr.integral_list) { total_integral += std::accumulate(inner_vec.begin(), inner_vec.end(), 0.0); }
    double sum_value = sr.u_tau_zeroth_order[sp.bl_index](i, j) + total_integral;
    double sum_time  = std::accumulate(sr.statistics[0].time_order.begin(), sr.statistics[0].time_order.end(), 0.0);
    std::cout << std::setw(10) << "sum:" << std::setw(30) << sum_value << std::setw(30) << sum_time << std::endl;

    std::cout << "---- results comparision ----" << std::endl;
    std::cout << "u_tau_max at tau=" << sp.tau_max << " for bl_index=" << sp.bl_index << ", subspace_index=" << sp.subspace_index << std::endl;

    // Table headers
    std::cout << std::setw(40) << std::left << "Description" << std::setw(30) << "Exact" << std::setw(30) << "Hyb" << std::endl;

    // u values
    std::cout << std::setw(40) << std::left << "u" << std::setw(30) << sr.u_interpolator_ref(sp.tau_max)[sp.bl_index](i, j) << std::setw(30)
              << sum_value * gp.Z_energy_shift_correction << std::endl;

    // u * Z_imp_correction values
    std::cout << std::setw(40) << std::left << "u*Z_imp_correction" << std::setw(30)
              << sr.u_interpolator_ref(sp.tau_max)[sp.bl_index](i, j) * sr.Z_imp_correction << std::setw(30)
              << sum_value * sr.Z_imp_correction * gp.Z_energy_shift_correction << std::endl;

    // u * Z_imp_correction * Z_bath * Z_bath_correction values
    std::cout << std::setw(40) << std::left << "u*Z_imp_correction*Z_bath_(correction)" << std::setw(30)
              << sr.u_interpolator_ref(sp.tau_max)[sp.bl_index](i, j) * sr.Z_bath * sr.Z_imp_correction * sr.Z_bath_correction << std::setw(30)
              << sum_value * sr.Z_bath * sr.Z_imp_correction * sr.Z_bath_correction * gp.Z_energy_shift_correction << std::endl;
  }
}

void ModeDebug::run() {
  if(rank == 0){
    std::cout << "### debug mode: start running ###" << std::endl;
  }
  validate_input();
  ModeBase::prepare_eval();
  if (gp.target == "propagator")
    evaluate_propagator();
  else if (gp.target == "greens_function")
    evaluate_greens_function();
  else
    std::cerr << "debug mode: invalid target" << std::endl;
}

void ModeDebug::evaluate_propagator() {
  sr.u_tau = sr.u_tau_ref;
  if (gp.do_regularization) {
    double exponent              = 0.0;
    std::tie(sr.u_tau, exponent) = regularize_propagator(sr.u_tau, mp, sp, sp.grid.size() - 1, gp.amplification_u);
    gp.exponent_u += exponent;
  }
  sr.u_interpolator      = interpolator_t<scalar_t>(sr.u_tau, sp.n_tot, sp.n_tau_linear, sp.order_Chebyshev, sp.interp_type, sp.grid);
  sr.u_tau_zeroth_order  = sr.u_interpolator(sp.tau_max - sp.tau_split) * sr.u_interpolator(sp.tau_split);
  sp.use_bare_propagator = false;
  ModeBase::evaluate();
}

void ModeDebug::evaluate_greens_function() {
  sr.u_tau = sr.u_tau_ref;
  ModeBase::evaluate_greens_function_bold();
}
