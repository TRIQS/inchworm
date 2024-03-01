#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "./mode.hpp"

using namespace inchworm;

void ModeInchworm::validate_input() { ModeBase::validate_input(); }

void ModeInchworm::run() {
  std::cout << "### Inchworm mode: start running ###" << std::endl;
  validate_input();
  evaluate_propagator();
}

void ModeInchworm::evaluate_propagator() {

  auto u_tau_zero = u_tau_t{{cp.beta, Fermion, cp.n_tau_inch}, mp.ad_imp.get_subspace_dims()};
  u_tau_zero()    = 0.;
  sr.u_tau        = u_tau_zero;
  for (auto &ubl : sr.u_tau) {
    for (int i = 0; i < ubl.target_shape()[0]; ++i) ubl[0](i, i) = 1;
  }
  double dtau = cp.beta / (cp.n_tau_inch - 1.);
  for (size_t i_tau = 1; i_tau < cp.n_tau_inch; i_tau++) {
    sp.tau_max             = sr.u_tau[0].mesh()(i_tau);
    sp.tau_split           = sr.u_tau[0].mesh()(i_tau - 1);
    std::cout << "sp.tau_max = " << sp.tau_max << std::endl;
    std::cout << "sp.tau_split = " << sp.tau_split << std::endl;
    sp.use_bare_propagator = (i_tau == 1);
    if (!sp.use_bare_propagator) { sr.u_interpolator = interpolator_t<scalar_t>(sr.u_tau, i_tau); }
    sr.u_tau_zeroth_order = sp.use_bare_propagator ? make_bare_u_frame(mp.ad_imp, sp.tau_max) :
                                                     sr.u_interpolator(sp.tau_max - sp.tau_split) * sr.u_interpolator(sp.tau_split);
    auto u_frame          = make_zero_frame(mp.ad_imp.get_subspace_dims());
    for (auto bl : range(mp.ad_imp.n_subspaces())) {
      for (auto i : range(mp.ad_imp.get_subspace_dim(bl))) {
        for (auto j : range(mp.ad_imp.get_subspace_dim(bl))) {
          std::cout << "bl = " << bl << ", i = " << i << ", j = " << j << std::endl;
          ModeBase::clear_tci_results();
          sp.bl_index = bl;
          sp.subspace_index = i * mp.ad_imp.get_subspace_dim(bl) + j;
          ModeBase::evaluate();
          u_frame[bl](i, j)  = sr.u_tau_zeroth_order[bl](i, j) + std::accumulate(sr.integral_list.begin(), sr.integral_list.end(), 0.0);
        } // end of j loop
      }   // end of j loop
    }     // end of bl loop

  set_frame(u_frame, sr.u_tau, i_tau);
  //print u_frame
  std::cout << "u_frame at tau = " << i_tau*dtau << std::endl;
      for (auto bl : range(mp.ad_imp.n_subspaces())) {
      for (auto i : range(mp.ad_imp.get_subspace_dim(bl))) {
        for (auto j : range(mp.ad_imp.get_subspace_dim(bl))) {
          std::cout << "u_frame[" << bl << "](" << i << "," << j << ") = " << u_frame[bl](i, j) << std::endl;
          std::cout << "sr.u_tau_ref[" << bl << "](" << i << "," << j << ") = " << sr.u_tau_ref[bl](i_tau*dtau)(i, j) << std::endl;
        } // end of j loop
      }   // end of j loop
    }     // end of bl loop

  } // end of i_tau loop
}

void ModeInchworm::evaluate_greens_function() {}