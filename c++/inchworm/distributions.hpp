#pragma once

#include "./types.hpp"
#include "./torus.hpp"

#include <triqs/mc_tools/random_generator.hpp>

#include <cmath>

namespace inchworm {

  // Truncated Cauchy Distribution on the interval [0,tmax)
  inline double pdf(double tau, double w, double tmax) {
    EXPECTS(tau >= 0 && tau <= tmax);
    return 1.0 / (w * (1 + tau * tau / (w * w)) * std::atan(tmax / w));
  }

  // Inverse of the cumulative distribution function
  // associated with pdf
  inline double icdf(double p, double w, double tmax) {
    EXPECTS(w > 0.0 && p >= 0.0);
    if (p >= 1.0) return tmax;
    double res = w * std::tan(p * std::atan(tmax / w));
    ASSERT(0.0 <= res);
    return res;
  }

  // Probability density function on the interval [0,tmax] combining
  // two probabily densities on [0,tmax/2) and [tmax/2,tmax]
  inline double double_pdf(double tau, double w1, double w2, double tmax) {
    EXPECTS(w1 > 0.0 && w2 > 0.0);
    EXPECTS(0 <= tau && tau <= tmax);
    // Piecewise pdf on [0,tmax/2) and [tmax/2,tmax]
    if (tau < tmax / 2) return 0.5 * pdf(tau, w1, tmax / 2.0);
    return 0.5 * pdf(tmax - tau, w2, tmax / 2.0);
  }

  // Inverse of the cumulative distribution function
  // associated with the double_pdf
  inline double double_icdf(double p, double w1, double w2, double tmax) {
    // Piecewise ipdf on [0,tmax/2) and [tmax/2,tmax]
    if (p < 0.5) return icdf(2.0 * p, w1, tmax / 2.0);
    return tmax - icdf(2.0 - 2.0 * p, w2, tmax / 2.0);
  }

  // -----------------------------------

  // Calculate the value of the joint pdf defined through split_times and the double_pdf function
  inline double get_prob(fop_t const &op, std::vector<double> const &split_times, double tmax) {
    double prob = 0.0;
    for (double tau_ref : split_times) {
      double diff = cyclic_difference(op.tau, tau_ref, tmax);
      prob += double_pdf(diff, op.left_width, op.right_width, tmax);
    }
    return prob / split_times.size();
  }

  // Draw a random time from the joint pdf defined through split_times the double_pdf function
  inline double get_close_time(auto &rng, fop_t const &op, std::vector<double> const &split_times, double tmax) {
    double tau_ref = split_times[rng(split_times.size())];
    return wrap(tau_ref + double_icdf(rng(), op.left_width, op.right_width, tmax), tmax);
  }

} // namespace inchworm
