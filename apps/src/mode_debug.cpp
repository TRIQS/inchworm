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
  ModeBase::print_summary();
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
  auto u_tau_ref_shifted = sr.u_tau_ref;
  if (gp.do_regularization) {
    std::tie(u_tau_ref_shifted, gp.exponent_u) = regularize_propagator(sr.u_tau_ref, mp, sp, sp.grid.size() - 1, gp.amplification_u);
  }
  auto u_interpolator_ref_shifted =
     interpolator_t<scalar_t>(u_tau_ref_shifted, sp.n_tot, sp.n_tau_linear, sp.order_Chebyshev, sp.interp_type, sp.grid);
  sr.u_tau_zeroth_order  = u_interpolator_ref_shifted(sp.tau_max - sp.tau_split) * u_interpolator_ref_shifted(sp.tau_split);
  long n_tot             = u_tau_ref_shifted[0].mesh().size();
  sr.u_interpolator      = interpolator_t<scalar_t>(u_tau_ref_shifted, sp.n_tot, sp.n_tau_linear, sp.order_Chebyshev, sp.interp_type, sp.grid);
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
  // treat the first and last point separately
  EXPECTS(sr.G_tau[0].mesh()[0] == 0.0);
  frame_t g_frame_0 = make_bare_g_frame(mp.ad_imp, sr.u_tau_ref, cp.gf_struct, 0.0, cp.beta) / Tr_Ubeta;
  set_frame(g_frame_0, sr.G_tau, 0);
  EXPECTS(sr.G_tau[0].mesh()[cp.n_tau_green - 1] == cp.beta);
  frame_t g_frame_beta = make_bare_g_frame(mp.ad_imp, sr.u_tau_ref, cp.gf_struct, cp.beta, cp.beta) / Tr_Ubeta;
  set_frame(g_frame_beta, sr.G_tau, cp.n_tau_green - 1);

  //set zero-th order
  for (size_t n = 1; n < cp.n_tau_green - 1; n++) {
    sp.tau_split            = sr.G_tau[0].mesh()[n];
    auto frame_zeroth_order = make_bare_g_frame(mp.ad_imp, sr.u_tau_ref, cp.gf_struct, sp.tau_split, cp.beta);
    for (int bl = 0; bl < mp.gf_block_shape.size(); bl++) {
      for (auto orb_d : range(mp.gf_block_shape[bl])) {
        for (auto orb_ddag : range(mp.gf_block_shape[bl])) {
          int orb_d_index                  = bl2_to_bl1_gf(bl, orb_d, mp.gf_block_shape);
          int orb_ddag_index               = bl2_to_bl1_gf(bl, orb_ddag, mp.gf_block_shape);
          sr.G_tau[bl][n](orb_d, orb_ddag) = (frame_zeroth_order[bl](orb_d, orb_ddag)) / Tr_Ubeta;
        }
      }
    }
  }
  if (gp.unsummed_tci == 0) {
    sp.tau_max = cp.beta;
    for (size_t n = 1; n < cp.n_tau_green - 1; n++) {
      std::cout << "evaluating tau[" << n << "] = " << sr.G_tau[0].mesh()[n] << std::endl;
      sp.tau_split = sr.G_tau[0].mesh()[n];
      // the most naive way is to loop over the spin-oribital index for d and ddag separately, However, since the Green's function is saved in block format, we only allow spin-orbital indices within the same block.
      for (int bl = 0; bl < mp.gf_block_shape.size(); bl++) {
        for (auto orb_d : range(mp.gf_block_shape[bl])) {
          for (auto orb_ddag : range(mp.gf_block_shape[bl])) {
            int orb_d_index    = bl2_to_bl1_gf(bl, orb_d, mp.gf_block_shape);
            int orb_ddag_index = bl2_to_bl1_gf(bl, orb_ddag, mp.gf_block_shape);
            sp.gf_index.clear();
            sp.gf_index.push_back(orb_d_index);
            sp.gf_index.push_back(orb_ddag_index);
            ModeBase::clear_tci_results();
            ModeBase::evaluate();
            double total_integral = 0.0;
            for (auto integral : sr.integral_list) total_integral += std::accumulate(integral.begin(), integral.end(), 0.0);
            double tci_result = total_integral / Tr_Ubeta;
            sr.G_tau[bl][n](orb_d, orb_ddag) += tci_result;
          }
        }
      }
    }
  } else if (gp.unsummed_tci == 1) {
    // include the time index into the tci
    std::vector<std::vector<double>> unsummed_input;
    std::vector<double> tau_list;
    for (size_t n = 1; n < cp.n_tau_green - 1; n++) { tau_list.push_back(sr.G_tau[0].mesh()[n]); }
    unsummed_input.push_back(tau_list);
    size_t dims_tau = tau_list.size();
    sp.tau_max      = cp.beta;
    sp.tau_split    = sr.G_tau[0].mesh()[1];
    for (int bl = 0; bl < mp.gf_block_shape.size(); bl++) {
      for (auto orb_d : range(mp.gf_block_shape[bl])) {
        for (auto orb_ddag : range(mp.gf_block_shape[bl])) {
          int orb_d_index    = bl2_to_bl1_gf(bl, orb_d, mp.gf_block_shape);
          int orb_ddag_index = bl2_to_bl1_gf(bl, orb_ddag, mp.gf_block_shape);
          sp.gf_index.clear();
          sp.gf_index.push_back(orb_d_index);
          sp.gf_index.push_back(orb_ddag_index);
          ModeBase::clear_tci_results();
          ModeBase::evaluate(unsummed_input);
          // save results
          std::vector<double> integrals(dims_tau, 0.0);
          for (auto integral_order : sr.integral_list) {
            for (size_t i = 0; i < integral_order.size(); i++) { integrals[i] += integral_order[i]; }
          }
          for (size_t n = 1; n < cp.n_tau_green - 1; n++) {
            auto tci_result = integrals[n - 1] / Tr_Ubeta;
            sr.G_tau[bl][n](orb_d, orb_ddag) += tci_result;
          }
        }
      }
    }
  } else if (gp.unsummed_tci == 2) {
    sp.gf_index.clear();
    sp.gf_index.push_back(0);
    sp.gf_index.push_back(0);
    //include the two orb indexs and tau index into tensor train
    std::vector<std::vector<double>> unsummed_input;
    size_t dims_orb = 0;
    for (int bl = 0; bl < mp.gf_block_shape.size(); bl++) { dims_orb += mp.gf_block_shape[bl] * mp.gf_block_shape[bl]; }
    std::vector<double> input(dims_orb, 0.0);
    std::iota(input.begin(), input.end(), 0);
    unsummed_input.push_back(input);
    std::vector<double> tau_list;
    for (size_t n = 1; n < cp.n_tau_green - 1; n++) { tau_list.push_back(sr.G_tau[0].mesh()[n]); }
    size_t dims_tau = tau_list.size();
    unsummed_input.push_back(tau_list);
    size_t dims_tot = dims_orb * dims_tau;
    sp.tau_max      = cp.beta;
    ModeBase::clear_tci_results();
    ModeBase::evaluate(unsummed_input);
    std::vector<double> integrals(dims_tot, 0.0);
    for (auto integral_order : sr.integral_list) {
      for (size_t i = 0; i < integral_order.size(); i++) { integrals[i] += integral_order[i]; }
    }
    for (size_t tau_index = 0; tau_index < dims_tau; tau_index++) {
      for (size_t orb_index = 0; orb_index < dims_orb; orb_index++) {
        size_t k                = tau_index + dims_tau * orb_index;
        auto [bl, i, j]         = bl1_to_bl3(orb_index, mp.gf_block_shape);
        size_t tau_index_actual = tau_index + 1;
        auto tci_result         = integrals[k] / Tr_Ubeta;
        sr.G_tau[bl][tau_index_actual](i, j) += tci_result;
      }
    }
  } else {
    throw std::runtime_error("evaluate Green's function: invalid unsummed_tci");
  }

  if (rank == 0) {
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
}
