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
  ModeBase::print_summary();


  auto [bl, i, j] = bl2_to_bl3(sp.bl_index, sp.subspace_index, mp.ad_imp.get_subspace_dims());
  std::cout << "#### debug results ####" << std::endl;
  std::cout << "---- order contribution ----" << std::endl;
  std::cout << std::left << std::setw(10) << "order" << std::setw(30) << "value" << std::setw(30) << "time(s)" << std::setw(30) << std::endl;
  std::cout << std::setw(10) << "0" << std::setw(30) << sr.u_tau_zeroth_order[sp.bl_index](i, j) << std::endl;
  for (int i = 0; i < sp.order_list.size(); i++) {
    std::cout << std::setw(10) << sp.order_list[i] << std::setw(30) << sr.integral_list[i][0] << std::setw(30) << (sr.statistics[0].time_order)[i]
              << std::endl;
  }
  double total_integral = 0.0;
  for (const auto &inner_vec : sr.integral_list) { total_integral += std::accumulate(inner_vec.begin(), inner_vec.end(), 0.0); }
  double sum_value = sr.u_tau_zeroth_order[sp.bl_index](i, j) + total_integral;
  double sum_time = std::accumulate(sr.statistics[0].time_order.begin(), sr.statistics[0].time_order.end(), 0.0);
    std::cout << std::setw(10) << "sum:" << std::setw(30) << sum_value << std::setw(30) << sum_time << std::endl;
  
  if(gp.target == "propagator"){
  std::cout << "---- results comparision ----" << std::endl;
  std::cout << "u_tau_max at tau=" << sp.tau_max << " for bl_index=" << sp.bl_index << ", subspace_index=" << sp.subspace_index << std::endl;

  // Table headers
  std::cout << std::setw(40) << std::left << "Description" 
            << std::setw(30) << "Exact" 
            << std::setw(30) << "Hyb" << std::endl;

  // u values
  std::cout << std::setw(40) << std::left << "u" 
            << std::setw(30) << sr.u_interpolator_ref(sp.tau_max)[sp.bl_index](i, j) 
            << std::setw(30) << sum_value * gp.Z_energy_shift_correction << std::endl;

  // u * Z_imp_correction values
  std::cout << std::setw(40) << std::left << "u*Z_imp_correction" 
            << std::setw(30) << sr.u_interpolator_ref(sp.tau_max)[sp.bl_index](i, j) * sr.Z_imp_correction 
            << std::setw(30) << sum_value * sr.Z_imp_correction * gp.Z_energy_shift_correction << std::endl;

  // u * Z_imp_correction * Z_bath * Z_bath_correction values
  std::cout << std::setw(40) << std::left << "u*Z_imp_correction*Z_bath_(correction)" 
            << std::setw(30) << sr.u_interpolator_ref(sp.tau_max)[sp.bl_index](i, j) * sr.Z_bath * sr.Z_imp_correction * sr.Z_bath_correction 
            << std::setw(30) << sum_value * sr.Z_bath * sr.Z_imp_correction * sr.Z_bath_correction * gp.Z_energy_shift_correction << std::endl;
  
  }
}

void ModeDebug::run() {
  std::cout << "### debug mode: start running ###" << std::endl;
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
  // for debug mode, the discrete bath is used
  sr.u_tau_zeroth_order = sr.u_interpolator_ref(sp.tau_max - sp.tau_split) * sr.u_interpolator_ref(sp.tau_split);
  long n_tot = sr.u_tau_ref[0].mesh().size();
  sr.u_interpolator =
     interpolator_t<scalar_t>(sr.u_tau_ref, n_tot, sp.n_tau_linear, sp.order_Chebyshev, interpolation_type::linear_Chebyshev, sp.grid);
  sp.use_bare_propagator = false;
  ModeBase::evaluate();
}

void ModeDebug::evaluate_greens_function() {
  // for debug mode, the discrete bath is used, so that we get u_interpolator from u_tau_ref
  sr.u_interpolator      = interpolator_t<scalar_t>(sr.u_tau_ref, sp.n_tot, sp.n_tau_linear, sp.order_Chebyshev, sp.interp_type, sp.grid);
  sp.use_bare_propagator = false;
  // Calculate Tr U(beta)
  scalar_t Tr_Ubeta = 0.0;
  for (int bl = 0; bl < sr.u_tau_ref.size(); bl++) Tr_Ubeta += trace(sr.u_tau_ref[bl][sp.n_tot - 1]);
  std::cout << "Z = Tr[U(beta)]: " << Tr_Ubeta << std::endl;

  sr.G_tau     = g_tau_t{{cp.beta, Fermion, cp.n_tau_green}, cp.gf_struct};
  int tot_dims = 0;
  for (auto subspace_dim : mp.gf_block_shape) { tot_dims += subspace_dim * subspace_dim; }
  for (size_t n = 0; n < cp.n_tau_green; n++) {
    std::cout << "evaluating tau[" << n << "] = " << sr.G_tau[0].mesh()[n] << std::endl;
    // treat the first and last point separately
    if (n == 0) {
      EXPECTS(sr.G_tau[0].mesh()[0] == 0.0);
      frame_t g_frame = make_bare_g_frame(mp.ad_imp, sr.u_tau_ref, cp.gf_struct, 0.0, cp.beta) / Tr_Ubeta;
      set_frame(g_frame, sr.G_tau, 0);
    } else if (n == cp.n_tau_green - 1) {
      EXPECTS(sr.G_tau[0].mesh()[cp.n_tau_green - 1] == cp.beta);
      frame_t g_frame = make_bare_g_frame(mp.ad_imp, sr.u_tau_ref, cp.gf_struct, cp.beta, cp.beta) / Tr_Ubeta;
      set_frame(g_frame, sr.G_tau, cp.n_tau_green - 1);
    } else { // other points
      sp.tau_split            = sr.G_tau[0].mesh()[n];
      sp.tau_max              = cp.beta;
      auto frame_zeroth_order = make_bare_g_frame(mp.ad_imp, sr.u_tau_ref, cp.gf_struct, sp.tau_split, cp.beta);
      for (int bl = 0; bl < mp.gf_block_shape.size(); bl++) {
        for (auto orb_d : range(mp.gf_block_shape[bl])) {
          for (auto orb_ddag : range(mp.gf_block_shape[bl])) {
            int orb_d_index    = bl2_to_bl1(orb_d, bl, mp.gf_block_shape);
            int orb_ddag_index = bl2_to_bl1(orb_ddag, bl, mp.gf_block_shape);
            sp.gf_index.clear();
            sp.gf_index.push_back(orb_d_index);
            sp.gf_index.push_back(orb_ddag_index);
            ModeBase::clear_tci_results();
            ModeBase::evaluate();
            double total_integral = 0.0;
            for (auto integral : sr.integral_list) total_integral += std::accumulate(integral.begin(), integral.end(), 0.0);
            double tci_result                = (frame_zeroth_order[bl](orb_d, orb_ddag) + total_integral) / Tr_Ubeta;
            sr.G_tau[bl][n](orb_d, orb_ddag) = tci_result;
          }
        }
      }
    }
  }
  // print the result
  for (size_t n = 0; n < cp.n_tau_green; n++) {
    for (int bl = 0; bl < mp.gf_block_shape.size(); bl++) {
      for (auto orb_d : range(mp.gf_block_shape[bl])) {
        for (auto orb_ddag : range(mp.gf_block_shape[bl])) {
          std::cout << "G[" << n << "][" << bl << "][" << orb_d << "][" << orb_ddag << "] = " << sr.G_tau[bl][n](orb_d, orb_ddag) << std::endl;
          std::cout << "G_ref[" << n << "][" << bl << "][" << orb_d << "][" << orb_ddag << "] = " << sr.G_tau_ref[bl][n](orb_d, orb_ddag)
                    << std::endl;
        }
      }
    }
  }

  auto file_name = gp.output_prefix + ".h5";
  h5::file file{file_name, 'w'};
  h5::group group{file};
  h5_save_params(this, group, "params");
  h5_save_gf(this, group, "gf", sr.G_tau);
  h5_save_gf(this, group, "gf_ref", sr.G_tau_ref);
}
