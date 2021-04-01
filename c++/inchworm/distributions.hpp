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
    EXPECTS(0.0 < std::exp(-tmax / w));
    // Explicitly treat small tmax / w to avoid instabilities
    if (tmax / w < 1e-8) return p * tmax;
    return -w * std::log(1.0 - p + p * std::exp(-tmax / w));
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

  static const std::pair<double, double> widths_d     = {2.0, 2.0};
  static const std::pair<double, double> widths_d_dag = {2.0, 2.0};

  // Calculate the value of the joint pdf defined through op_list and the double_pdf function
  inline double get_prob(double tau, std::vector<fop_t> const &op_list, double tmax) {
    double prob = 0.0;
    for (auto const &op : op_list) {
      double diff   = cyclic_difference(tau, op.tau, tmax);
      auto [w1, w2] = op.dag ? widths_d : widths_d_dag;
      prob += double_pdf(diff, w1, w2, tmax);
    }
    return prob / op_list.size();
  }

  // Draw a random time from the joint pdf defined through op_list the double_pdf function
  inline double get_close_time(triqs::mc_tools::random_generator &rng, std::vector<fop_t> const &op_list, double tmax) {
    auto const &op = op_list[rng(op_list.size())];
    auto [w1, w2]  = op.dag ? widths_d : widths_d_dag;
    return wrap(op.tau + double_icdf(rng(), w1, w2, tmax), tmax);
  }

  // Draw a random time from the joint pdf defined through op_list the double_pdf function
  // and return both the time and its proposition probability
  inline std::pair<double, double> get_close_time_and_prob(triqs::mc_tools::random_generator &rng, std::vector<fop_t> const &op_list, double tmax) {
    double tau = get_close_time(rng, op_list, tmax);
    return {tau, get_prob(tau, op_list, tmax)};
  }

} // namespace inchworm
