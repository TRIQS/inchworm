#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "./mode.hpp"

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

void ModeDebug::run() {
  std::cout << "### debug mode: start running ###" << std::endl;
  validate_input();
  if (gp.target == "propagator")
    evaluate_propagator();
  else if (gp.target == "greens_function")
    evaluate_greens_function();
  else
    std::cerr << "debug mode: invalid target" << std::endl;
}

void ModeDebug::evaluate_propagator() {
  // for debug mode, the discrete bath is used
  sr.u_tau_zeroth_order = sr.u_tau_zeroth_order_ref;
  // sr.u_interpolator      = interpolator_t<scalar_t>(sr.u_tau_ref, sr.u_tau_ref[0].mesh().size(), 0, 0, interpolation_type::cspline);
  long n_tot = sr.u_tau_ref[0].mesh().size();
  std::cout << "n_tot: " << n_tot << ", n_tau: " << sp.n_tau_linear << ", order: " << sp.order_Chebyshev << std::endl;
  sr.u_interpolator      = interpolator_t<scalar_t>(sr.u_tau_ref, n_tot, sp.n_tau_linear, sp.order_Chebyshev, interpolation_type::linear_Chebyshev);
  std::cout<<"grid u_tau_ref:"<<std::endl;
  for(auto tau: sr.u_tau_ref[0].mesh()){
    std::cout<<tau<<std::endl;
  }
  // auto [grid_linear, grid] = generate_linear_Chebyshev_grid(0, cp.beta, sp.n_tau_linear, sp.order_Chebyshev);
  // std::cout<<"grid: ";
  // print_vector(grid);
  sp.use_bare_propagator = false;
  ModeBase::evaluate();
}

void ModeDebug::evaluate_greens_function() {
  // for debug mode, the discrete bath is used
  sr.u_tau_zeroth_order  = sr.u_tau_zeroth_order_ref;
  sr.u_interpolator      = interpolator_t<scalar_t>(sr.u_tau_ref, sr.u_tau_ref[0].mesh().size(), 0, 0, interpolation_type::cspline);
  sp.use_bare_propagator = false;
  sp.eval_type           = 1;
  sp.n_skip              = 0;

  sr.G_tau = g_tau_t{{cp.beta, Fermion, cp.n_tau_green}, cp.gf_struct};
  // Calculate Tr U(beta)
  scalar_t Tr_Ubeta = 0.0;
  for (int bl = 0; bl < sr.u_tau_ref.size(); bl++) Tr_Ubeta += trace(sr.u_tau_ref[bl](cp.beta));

  // Treat n == 0 and n == n_tau -1 seperately
  frame_t g_frame_n0 = make_bare_g_frame(mp.ad_imp, sr.u_tau_ref, cp.gf_struct, 0.0, cp.beta) / Tr_Ubeta;
  frame_t g_frame_nB = make_bare_g_frame(mp.ad_imp, sr.u_tau_ref, cp.gf_struct, cp.beta, cp.beta) / Tr_Ubeta;
  set_frame(g_frame_n0, sr.G_tau, 0);
  set_frame(g_frame_nB, sr.G_tau, cp.n_tau_green - 1);
  std::cout << "sr.G_tau[0][0](0,0): " << sr.G_tau[0][0](0, 0) << std::endl;
  std::cout << "sr.G_tau[0][beta](0,0): " << sr.G_tau[0][cp.n_tau_green - 1](0, 0) << std::endl;
  for (size_t n = 0; n < cp.n_tau_green; n++) {
    std::cout << "n: " << n << std::endl;
    std::cout << "tau[n]: " << sr.G_tau[0].mesh()[n] << std::endl;
    std::cout << "sr.G_tau[0][n](0,0): " << sr.G_tau[0][n](0, 0) << std::endl;
  }
  sp.tau_split = 0.5 * cp.beta;
  sp.gf_index.push_back(0);
  sp.gf_index.push_back(0);
  ModeBase::evaluate();
  auto frame_zeroth_order = make_bare_g_frame(mp.ad_imp, sr.u_tau_ref, cp.gf_struct, sp.tau_split, cp.beta);
  double tci_result       = (frame_zeroth_order[0](0, 0) + std::accumulate(sr.integral_list.begin(), sr.integral_list.end(), 0.0)) / Tr_Ubeta;
  std::cout << "G(tau_split): " << sr.G_tau_ref[0][1](0, 0) << std::endl;
  std::cout << "tci_result: " << tci_result << std::endl;
}