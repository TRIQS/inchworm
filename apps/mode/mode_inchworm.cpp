#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "./mode.hpp"
#include "../save.hpp"

using namespace inchworm;

void ModeInchworm::validate_input() { ModeBase::validate_input(); }

void ModeInchworm::run() {
  std::cout << "### Inchworm mode: start running ###" << std::endl;
  validate_input();
  evaluate_propagator();
}

void ModeInchworm::evaluate_propagator() {
  // direct sum tests; the first layer is for different order, the second layer is for different index in TT, the third layer is for different indices for a given index in TT
  std::vector<std::vector<std::vector<double>>> all_input{};
  std::vector<std::vector<double>> all_weight{};
  if (gp.exact_sum) {
    auto [input_order1, weight_order1] = generate_combination_at_fixed_order(1, tp.v_value, tp.v_weight);
    auto [input_order2, weight_order2] = generate_combination_at_fixed_order(2, tp.v_value, tp.v_weight);
    all_input.push_back(input_order1);
    all_input.push_back(input_order2);
    all_weight.push_back(weight_order1);
    all_weight.push_back(weight_order2);
  }

  std::vector<std::vector<double>> unsummed_input{};
  if (gp.unsummed_tci) {
    int total_dims = 0;
    for (auto bl : range(mp.ad_imp.n_subspaces())) { total_dims += mp.ad_imp.get_subspace_dim(bl); }
    std::vector<double> input(total_dims, 0.0);
    std::iota(input.begin(), input.end(), 0.0);
    unsummed_input.push_back(std::vector<double>{input});
  }

  auto u_tau_zero = u_tau_t{{cp.beta, Fermion, sp.n_tot}, mp.ad_imp.get_subspace_dims()};
  u_tau_zero()    = 0.;
  sr.u_tau        = u_tau_zero;
  for (auto &ubl : sr.u_tau) {
    for (int i = 0; i < ubl.target_shape()[0]; ++i) ubl[0](i, i) = 1;
  }
  // inchworm grid points are in sp.grid_linear of length sp.n_tau_linear
  // 0 [0], 1 [dtau], ... , sp.n_tau_linear-1 [beta]
  // all the evaluated points are in sp.grid of length sp.n_tot
  for (size_t i_tau = 1; i_tau < sp.n_tau_linear; i_tau++) {
    sp.use_bare_propagator     = (i_tau == 1);
    long i_grid_tau_split      = i_tau + (i_tau - 1) * (sp.order_Chebyshev + 1) - 1;
    long i_grid_tau_next_split = i_tau + 1 + (i_tau) * (sp.order_Chebyshev + 1) - 1;
    EXPECTS(sp.grid[i_grid_tau_split] == sp.grid_linear[i_tau - 1]);
    EXPECTS(sp.grid[i_grid_tau_next_split] == sp.grid_linear[i_tau]);
    if (!sp.use_bare_propagator) {
      std::vector<double> grid_tau_split(sp.grid.begin(), sp.grid.begin() + i_grid_tau_split + 1);
      sr.u_interpolator =
         interpolator_t<scalar_t>(sr.u_tau, i_grid_tau_split + 1, i_tau, sp.order_Chebyshev, interpolation_type::linear_Chebyshev, grid_tau_split);
    }
    sp.tau_split = sp.grid_linear[i_tau - 1];
    // evaluate points from sp.grid[i_tau+(i_tau-1)*(order_Chebyshev+1)] to sp.grid[i_tau+1+(i_tau)*(order_Chebyshev+1)-1]
    for (size_t i_Chebyshev_tau = 0; i_Chebyshev_tau < sp.order_Chebyshev + 2; i_Chebyshev_tau++) {
      sp.tau_max = sp.grid[i_grid_tau_split + 1 + i_Chebyshev_tau];
      std::cout << "sp.tau_max = " << sp.tau_max << std::endl;
      std::cout << "sp.tau_split = " << sp.tau_split << std::endl;
      sr.u_tau_zeroth_order = sp.use_bare_propagator ? make_bare_u_frame(mp.ad_imp, sp.tau_max) :
                                                       sr.u_interpolator(sp.tau_max - sp.tau_split) * sr.u_interpolator(sp.tau_split);
      auto u_frame          = make_zero_frame(mp.ad_imp.get_subspace_dims());

      if (gp.unsummed_tci) {
        ModeBase::clear_tci_results();
        ModeBase::evaluate(unsummed_input, all_input, all_weight);
        int total_dims = 0;
        for (auto bl : range(mp.ad_imp.n_subspaces())) { total_dims += mp.ad_imp.get_subspace_dim(bl); }
        std::vector<double> integrals(total_dims, 0.0);
        for (auto integral_order : sr.integral_list){
          for(size_t i = 0; i < integral_order.size(); i++){
            integrals[i] += integral_order[i];
          }
        }
      for(size_t k = 0; k < integrals.size(); k++){
        auto[bl,i,j]= bl1_to_bl3(k, mp.ad_imp.get_subspace_dims());
        u_frame[bl](i,j) = sr.u_tau_zeroth_order[bl](i,j) + integrals[k];
      }
      } else {
        for (auto bl : range(mp.ad_imp.n_subspaces())) {
          for (auto i : range(mp.ad_imp.get_subspace_dim(bl))) {
            for (auto j : range(mp.ad_imp.get_subspace_dim(bl))) {
              std::cout << "bl = " << bl << ", i = " << i << ", j = " << j << std::endl;
              ModeBase::clear_tci_results();
              sp.bl_index       = bl;
              sp.subspace_index = i * mp.ad_imp.get_subspace_dim(bl) + j;
              ModeBase::evaluate(unsummed_input, all_input, all_weight);
              // ModeBase::evaluate();
              double total_integral = 0.;
              for (auto integral : sr.integral_list) { total_integral += std::accumulate(integral.begin(), integral.end(), 0.0); }
              u_frame[bl](i, j) = sr.u_tau_zeroth_order[bl](i, j) + total_integral;
            } // end of j loop
          } // end of j loop
        } // end of bl loop
      }
      set_frame(u_frame, sr.u_tau, i_grid_tau_split + 1 + i_Chebyshev_tau);
      //print u_frame
      std::cout << "u_frame at tau = " << sp.tau_max << std::endl;
      for (auto bl : range(mp.ad_imp.n_subspaces())) {
        for (auto i : range(mp.ad_imp.get_subspace_dim(bl))) {
          for (auto j : range(mp.ad_imp.get_subspace_dim(bl))) {
            std::cout << "u_frame[" << bl << "](" << i << "," << j << ") = " << u_frame[bl](i, j) << std::endl;
            if (gp.model_type == 0) {
              std::cout << "sr.u_tau_ref[" << bl << "](" << i << "," << j << ") = " << sr.u_interpolator_ref(sp.tau_max)[bl](i, j) << std::endl;
            }
          } // end of j loop
        } // end of j loop
      } // end of bl loop
    }
  } // end of i_tau loop
  std::cout << "---- partition function ----" << std::endl;
  double partition_function     = 0.;
  double partition_function_ref = 0.;
  for (auto bl : range(mp.ad_imp.n_subspaces())) {
    for (auto i : range(mp.ad_imp.get_subspace_dim(bl))) {
      partition_function += sr.u_tau[bl][sp.n_tot - 1](i, i);
      if (gp.model_type == 0) { partition_function_ref += sr.u_tau_ref[bl][sp.n_tot - 1](i, i); }
    }
  }
  std::cout << "partition_function = " << partition_function << std::endl;
  std::cout << "partition_function_ref = " << partition_function_ref << std::endl;

  sr.u_interpolator =
     interpolator_t<scalar_t>(sr.u_tau, sp.n_tot, sp.n_tau_linear, sp.order_Chebyshev, interpolation_type::linear_Chebyshev, sp.grid);

  auto file_name = gp.output_prefix + ".h5";
  h5::file file{file_name, 'w'};
  h5::group group{file};
  h5_save_params(this, group, "params");
  h5_save_propagator(this, group, "propagator");
  h5_save_cheb_coeff(this, group, "cheb_coeff");
}

void ModeInchworm::evaluate_greens_function() {}