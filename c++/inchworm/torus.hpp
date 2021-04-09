#pragma once

#include "./types.hpp"
#include "./config.hpp"

namespace inchworm {

  // Find the standard representative of tau on the torus R + Z * tau_max
  inline double wrap(double tau, double tau_max) { return (tau > 0) ? std::fmod(tau, tau_max) : tau_max - std::fmod(-tau, tau_max); }

  // Return the difference of the times on the torus R + Z * tau_max
  inline double cyclic_difference(double tau1, double tau2, double tau_max) { return wrap(tau1 - tau2, tau_max); }

  // Return the shortest absolute distance of the times on the torus R + Z * tau_max
  inline double cyclic_distance(double tau1, double tau2, double tau_max) {
    auto wrap_abs = wrap(std::abs(tau1 - tau2), tau_max);
    return std::min(wrap_abs, tau_max - wrap_abs);
  }

  // Return the cyclic_difference between the creation/annihilation operator
  // and the closest annihilation/creation operator of the same block of the config
  inline double cyclic_difference(auto const &op, config_t const &config, double tau_max) {
    double min_dist = tau_max;
    double res      = tau_max;
    auto check_time_diff = [&](double tau) {
      auto diff = cyclic_difference(op.tau, tau, tau_max);
      auto dist = cyclic_distance(op.tau, tau, tau_max);
      if (dist < min_dist) {
        min_dist = dist;
        res      = diff;
      }
    };

    for (double tau_ref: config.split_times) check_time_diff(tau_ref);

    return res;
  }

} // namespace inchworm
