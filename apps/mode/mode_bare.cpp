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

void ModeBare::evaluate_greens_function() {}