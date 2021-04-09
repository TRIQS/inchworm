#pragma once

#include "./types.hpp"
#include "./torus.hpp"

#include <triqs/mc_tools/random_generator.hpp>

namespace inchworm {

  // Exponential probability density truncated to the interval [0,tmax)
  inline double pdf(double tau, double w, double tmax) { return std::exp(-tau / w) / (1.0 - std::exp(-tmax / w)) / w; }

  // Inverse of the cumulative distribution function associated
  // with the truncated exponential probability density
  inline double icdf(double p, double w, double tmax) {
    // Explicitly treat small tmax / w to avoid instabilities
    if (tmax / w < 1e-8) return p * tmax;
    double exp_val = std::exp(-tmax / w);
    if ((1.0 - p) * 1e+3 <= exp_val) return tmax;
    return -w * std::log(1.0 - p + p * exp_val);
  }

  // Probability density function on the interval [0,tmax] combining
  // two exponential probabily densities on [0,tmax/2) and [tmax/2,tmax]
  inline double double_pdf(double tau, double w1, double w2, double tmax) {
    EXPECTS(0 <= tau && tau <= tmax);

    // Piecewise pdf on [0,tmax/2) and [tmax/2,tmax]
    if (tau < tmax / 2) return 0.5 * pdf(tau, w1, tmax / 2.0);
    return 0.5 * pdf(tmax - tau, w2, tmax / 2.0);
  }

  // Inverse of the cumulative distribution function associated
  // with the truncated exponential probability density
  inline double double_icdf(double p, double w1, double w2, double tmax) {
    // Piecewise ipdf on [0,tmax/2) and [tmax/2,tmax]
    if (p < 0.5) return icdf(2.0 * p, w1, tmax / 2.0);
    return tmax - icdf(2.0 - 2.0 * p, w2, tmax / 2.0);
  }

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
  inline double get_close_time(triqs::mc_tools::random_generator &rng, fop_t const &op, std::vector<double> const &split_times, double tmax) {
    double tau_ref = split_times[rng(split_times.size())];
    return wrap(tau_ref + double_icdf(rng(), op.left_width, op.right_width, tmax), tmax);
  }

  // Draw a random time from the joint pdf defined through split_times the double_pdf function
  // and return both the time and its proposition probability
  inline std::pair<double, double> get_close_time_and_prob(triqs::mc_tools::random_generator &rng, fop_t op, std::vector<double> const &split_times,
                                                           double tmax) {
    op.tau = get_close_time(rng, op, split_times, tmax);
    return {op.tau, get_prob(op, split_times, tmax)};
  }

} // namespace inchworm
