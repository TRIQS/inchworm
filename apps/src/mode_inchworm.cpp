#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "mode.hpp"
#include "save.hpp"

using namespace inchworm;

void ModeInchworm::validate_input() { ModeBase::validate_input(); }

void ModeInchworm::print_summary() {
  if (sp.debug <= 0) return;
  ModeBase::print_params();
}

void ModeInchworm::run() {
  NVTX_RANGE("inchworm run", 3);
  if (rank == 0) { std::cout << "### Inchworm mode: start running ###" << std::endl; }
  validate_input();
  ModeBase::prepare_eval();
  if (gp.target == "propagator") {
    evaluate_propagator();
  } else if (gp.target == "greens_function") {
    gp.target = "propagator";
    evaluate_propagator();
    gp.target = "greens_function";
    evaluate_greens_function();
  } else if (gp.target == "greens_function_restart") {
    read_propagator();
    std::cout << "read propagator from file" << std::endl;
    gp.target = "greens_function";
    evaluate_greens_function();
  } else {
    std::cerr << "invalid target" << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

void ModeInchworm::evaluate_propagator() {

  std::vector<std::vector<double>> unsummed_input{};
  if (gp.unsummed_tci == 1 || gp.unsummed_tci == 2) {
    int total_dims = 0;
    for (auto bl : range(mp.ad_imp.n_subspaces())) { total_dims += mp.ad_imp.get_subspace_dim(bl) * mp.ad_imp.get_subspace_dim(bl); }
    std::vector<double> input(total_dims, 0.0);
    std::iota(input.begin(), input.end(), 0.0);
    unsummed_input.push_back(input);
  }

  auto u_tau_zero = u_tau_t{{cp.beta, Fermion, sp.n_tot}, mp.ad_imp.get_subspace_dims()};
  u_tau_zero()    = 0.;
  sr.u_tau        = u_tau_zero;
  for (auto &ubl : sr.u_tau) {
    for (int i = 0; i < ubl.target_shape()[0]; ++i) ubl[0](i, i) = 1;
  }

  if (sp.inch_start_index < 1 || sp.inch_start_index >= sp.n_tau_linear) {
    std::cerr << "Error: inch_start_index should be between 1 and n_tau_linear-1" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (sp.inch_start_index != 1 && gp.model_type != 0) {
    std::cerr << "Error: inch_start_index should be 1 for model_type != 0" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (sp.inch_end_index != (sp.n_tau_linear - 1) && gp.model_type != 0) {
    std::cerr << "Error: inch_end_index should be n_tau_linear-1 for model_type != 0" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (sp.inch_end_index < sp.inch_start_index || sp.inch_end_index >= sp.n_tau_linear) {
    std::cerr << "Error: inch_end_index should be between inch_start_index and n_tau_linear-1" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (sp.inch_start_index != 1) {
    for (int i = 0; i < (sp.inch_start_index - 1) * (sp.order_Chebyshev + 1) + sp.inch_start_index; i++) {
      for (int bl = 0; bl < sr.u_tau_ref.size(); bl++) sr.u_tau[bl][i] = sr.u_tau_ref[bl][i];
    }
  }

  // inchworm grid points are in sp.grid_linear of length sp.n_tau_linear
  // 0 [0], 1 [dtau], ... , sp.n_tau_linear-1 [beta]
  // all the evaluated points are in sp.grid of length sp.n_tot
  for (size_t i_tau = sp.inch_start_index; i_tau < sp.n_tau_linear; i_tau++) {
    if (i_tau > sp.inch_end_index) { break; }
    bool is_first_interval     = (i_tau == 1);
    sp.use_bare_propagator     = (i_tau == 1);
    long i_grid_tau_split      = i_tau + (i_tau - 1) * (sp.order_Chebyshev + 1) - 1;
    long i_grid_tau_next_split = i_tau + 1 + (i_tau) * (sp.order_Chebyshev + 1) - 1;
    EXPECTS(sp.grid[i_grid_tau_split] == sp.grid_linear[i_tau - 1]);
    EXPECTS(sp.grid[i_grid_tau_next_split] == sp.grid_linear[i_tau]);
    if (!sp.use_bare_propagator) {
      if (gp.do_regularization) {
        double exponent              = 0.0;
        std::tie(sr.u_tau, exponent) = regularize_propagator(sr.u_tau, mp, sp, i_grid_tau_split, gp.amplification_u);
        gp.exponent_u += exponent;
      }
      std::vector<double> grid_tau_split(sp.grid.begin(), sp.grid.begin() + i_grid_tau_split + 1);
      sr.u_interpolator =
         interpolator_t<scalar_t>(sr.u_tau, i_grid_tau_split + 1, i_tau, sp.order_Chebyshev, interpolation_type::linear_Chebyshev, grid_tau_split);
    }
    sp.tau_split = sp.grid_linear[i_tau - 1];
    // evaluate points from sp.grid[i_tau+(i_tau-1)*(order_Chebyshev+1)] to sp.grid[i_tau+1+(i_tau)*(order_Chebyshev+1)-1]

    if (gp.unsummed_tci != 2) {
      for (size_t i_Chebyshev_tau = 0; i_Chebyshev_tau < sp.order_Chebyshev + 2; i_Chebyshev_tau++) {
        sp.tau_max            = sp.grid[i_grid_tau_split + 1 + i_Chebyshev_tau];
        sr.u_tau_zeroth_order = sp.use_bare_propagator ? make_bare_u_frame(mp.ad_imp, sp.tau_max, gp.exponent_u) :
                                                         sr.u_interpolator(sp.tau_max - sp.tau_split) * sr.u_interpolator(sp.tau_split);
        auto u_frame          = make_zero_frame(mp.ad_imp.get_subspace_dims());

        if (gp.unsummed_tci == 1) {
          ModeBase::clear_tci_results();
          ModeBase::evaluate(unsummed_input, is_first_interval, i_tau - 1);
          int total_dims = 0;
          for (auto bl : range(mp.ad_imp.n_subspaces())) { total_dims += mp.ad_imp.get_subspace_dim(bl) * mp.ad_imp.get_subspace_dim(bl); }
          std::vector<double> integrals(total_dims, 0.0);
          for (auto integral_order : sr.integral_list) {
            for (size_t i = 0; i < integral_order.size(); i++) { integrals[i] += integral_order[i]; }
          }
          for (size_t k = 0; k < integrals.size(); k++) {
            auto [bl, i, j]   = bl1_to_bl3(k, mp.ad_imp.get_subspace_dims());
            u_frame[bl](i, j) = sr.u_tau_zeroth_order[bl](i, j) + integrals[k];
          }
        } else {
          for (auto bl : range(mp.ad_imp.n_subspaces())) {
            for (auto i : range(mp.ad_imp.get_subspace_dim(bl))) {
              for (auto j : range(mp.ad_imp.get_subspace_dim(bl))) {
                sp.bl_index       = bl;
                sp.subspace_index = i * mp.ad_imp.get_subspace_dim(bl) + j;
                ModeBase::evaluate(unsummed_input, is_first_interval, i_tau - 1);
                double total_integral = 0.;
                for (auto integral : sr.integral_list) { total_integral += std::accumulate(integral.begin(), integral.end(), 0.0); }
                u_frame[bl](i, j) = sr.u_tau_zeroth_order[bl](i, j) + total_integral;
              } // end of j loop
            } // end of j loop
          } // end of bl loop
        }
        set_frame(u_frame, sr.u_tau, i_grid_tau_split + 1 + i_Chebyshev_tau);
      }
    } else {
      // unsummed_tci == 2; we add the tmax index into the tensor train
      int dims_tau = sp.order_Chebyshev + 2;
      sp.tau_max   = sp.grid[i_grid_tau_split + 2 + sp.order_Chebyshev];
      std::vector<double> input(dims_tau, 0.0);
      for (size_t i_Chebyshev_tau = 0; i_Chebyshev_tau < sp.order_Chebyshev + 2; i_Chebyshev_tau++) {
        input[i_Chebyshev_tau] = sp.grid[i_grid_tau_split + 1 + i_Chebyshev_tau];
      }
      if (unsummed_input.size() == 2) { unsummed_input.pop_back(); }
      unsummed_input.push_back(input);
      assert(unsummed_input.size() == 2);
      ModeBase::clear_tci_results();
      ModeBase::evaluate(unsummed_input, is_first_interval, i_tau - 1);
      int dims_orb = 0;
      for (auto bl : range(mp.ad_imp.n_subspaces())) { dims_orb += mp.ad_imp.get_subspace_dim(bl) * mp.ad_imp.get_subspace_dim(bl); }
      std::vector<double> integrals(dims_tau * dims_orb, 0.0);
      for (auto integral_order : sr.integral_list) {
        for (size_t i = 0; i < integral_order.size(); i++) { integrals[i] += integral_order[i]; }
      }
      for (size_t tau_index = 0; tau_index < dims_tau; tau_index++) {
        auto u_frame = make_zero_frame(mp.ad_imp.get_subspace_dims());
        for (size_t orb_index = 0; orb_index < dims_orb; orb_index++) {
          size_t k              = tau_index + dims_tau * orb_index;
          auto [bl, i, j]       = bl1_to_bl3(orb_index, mp.ad_imp.get_subspace_dims());
          double tau_max        = sp.grid[i_grid_tau_split + 1 + tau_index];
          sr.u_tau_zeroth_order = sp.use_bare_propagator ? make_bare_u_frame(mp.ad_imp, tau_max, gp.exponent_u) :
                                                           sr.u_interpolator(tau_max - sp.tau_split) * sr.u_interpolator(sp.tau_split);
          u_frame[bl](i, j)     = sr.u_tau_zeroth_order[bl](i, j) + integrals[k];
        }
        set_frame(u_frame, sr.u_tau, i_grid_tau_split + 1 + tau_index);
      }
    }
    if (sp.debug > 0 && rank == 0) {
      std::cout << "---- inchworm ----" << std::endl;
      std::cout << "gp.exponent_u: " << gp.exponent_u << std::endl;
      std::cout << "tau_split = " << sp.tau_split << std::endl;
      std::cout << "tau_max = " << std::endl;
      if (gp.unsummed_tci == 2) {
        for (auto tau : unsummed_input[1]) { std::cout << tau << " "; }
      } else {
        for (size_t i_Chebyshev_tau = 0; i_Chebyshev_tau < sp.order_Chebyshev + 2; i_Chebyshev_tau++) {
          double tau_max = sp.grid[i_grid_tau_split + 1 + i_Chebyshev_tau];
          std::cout << tau_max << std::endl;
        }
      }
      std::cout << std::endl;
    }
  } // end of i_tau loop
  // set sr.u_tau be sr.u_tau_ref for i_tau > sp.inch_end_index
  if (sp.inch_end_index < sp.n_tau_linear - 1) {
    if (gp.model_type == 0) {
      for (int i = (sp.inch_end_index) * (sp.order_Chebyshev + 1) + sp.inch_end_index + 1; i < sp.n_tot; i++) {
        for (int bl = 0; bl < sr.u_tau_ref.size(); bl++) sr.u_tau[bl][i] = sr.u_tau_ref[bl][i] * std::exp(-gp.exponent_u * sp.grid[i]);
      }
    } else {
      throw std::runtime_error("inchworm mode does not support model_type != 0 while sp.inch_end_index < n_tau_linear - 1");
    }
  }
  if (gp.do_regularization) {
    double exponent              = 0.0;
    std::tie(sr.u_tau, exponent) = regularize_propagator(sr.u_tau, mp, sp, sp.n_tot - 1, gp.amplification_u);
    gp.exponent_u += exponent;
  }
  gp.Z_energy_shift_correction = std::exp(gp.exponent_u * cp.beta);
  if (sp.debug > 0 && rank == 0) {
    std::cout << "---- partition function ----" << std::endl;
    std::cout << "gp.Z_energy_shift_correction = " << gp.Z_energy_shift_correction << std::endl;
    double partition_function     = 0.;
    double partition_function_ref = 0.;
    for (auto bl : range(mp.ad_imp.n_subspaces())) {
      for (auto i : range(mp.ad_imp.get_subspace_dim(bl))) {
        partition_function += sr.u_tau[bl][sp.n_tot - 1](i, i);
        if (gp.model_type == 0) { partition_function_ref += sr.u_tau_ref[bl][sp.n_tot - 1](i, i); }
      }
    }
    std::cout << "partition_function*gp.Z_energy_shift_correction  = " << partition_function * gp.Z_energy_shift_correction << std::endl;
    std::cout << "partition_function_ref = " << partition_function_ref << std::endl;
  }

  sr.u_interpolator =
     interpolator_t<scalar_t>(sr.u_tau, sp.n_tot, sp.n_tau_linear, sp.order_Chebyshev, interpolation_type::linear_Chebyshev, sp.grid);
  if (rank == 0) {
    std::cout << "---- saving results ----" << std::endl;
    auto file_name = gp.output_prefix + ".h5";
    h5::file file{file_name, 'w'};
    h5::group group{file};
    h5_save_params(this, group, "params");
    h5_save_propagator(this, group, "propagator");
    h5_save_cheb_coeff(this, group, "cheb_coeff");
    if(gp.model_type == 0){
      h5_save_propagator_ref(this, group, "propagator_ref");
      h5_save_cheb_coeff_ref(this, group, "cheb_coeff_ref");
    }
    h5_save_statistics(this, group, "statistics");
  }
}

void ModeInchworm::evaluate_greens_function() { ModeBase::evaluate_greens_function_bold(); }