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
    EXPECTS(w > 0.0 && p >= 0.0);
    // Explicitly treat small tmax / w to avoid instabilities
    if (p >= 1.0) return tmax;
    if (tmax / w < 1e-8) return p * tmax;
    double res = std::min(-w * std::log(1.0 - p + p * std::exp(-tmax / w)), tmax);
    ASSERT(0.0 <= res);
    return res;
  }

  // Probability density function on the interval [0,tmax] combining
  // two exponential probabily densities on [0,tmax/2) and [tmax/2,tmax]
  inline double double_pdf(double tau, double w1, double w2, double tmax) {
    EXPECTS(w1 > 0.0 && w2 > 0.0);
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
  inline double get_prob(fop_t const &op, std::vector<fop_t> const &op_list, double tmax) {
    double prob = 0.0;
    for (auto const &op_ref : op_list) {
      double diff = cyclic_difference(op.tau, op_ref.tau, tmax);
      prob += double_pdf(diff, op.left_width, op.right_width, tmax);
    }
    return prob / op_list.size();
  }

  // Calculate the probability to draw op.tau larger than tref
  inline double get_prob_larger_time(fop_t const &op, double tref, double range, double period) {
    double diff = cyclic_difference(op.tau, tref, period);
    if (diff <= range) return pdf(diff, op.left_width, range);
    return 0.0;
  }

  // Calculate the probability to draw op.tau smaller than tref
  inline double get_prob_smaller_time(fop_t const &op, double tref, double range, double period) {
    double diff = cyclic_difference(tref, op.tau, period);
    if (diff <= range) return pdf(diff, op.right_width, range);
    return 0.0;
  }

  // Calculate the value of the joint pdf defined through op_list and the double_pdf function
  inline double get_prob_smaller_and_larger_time(fop_t const &op1, fop_t const &op2, std::vector<double> const &split_times, double tmax) {
    double prob = 0.0;
    for (auto const &tref : split_times) {
      prob += get_prob_smaller_time(op1, tref, tmax / 2.0, tmax) * get_prob_larger_time(op2, tref, tmax / 2.0, tmax);
      prob += get_prob_larger_time(op1, tref, tmax / 2.0, tmax) * get_prob_smaller_time(op2, tref, tmax / 2.0, tmax);
    }
    return prob / split_times.size() / 2.0;
  }

  // Calculate the value of the joint pdf defined through op_list and the double_pdf function
  inline double get_prob_smaller_and_larger_time(fop_t const &op1, fop_t const &op2, double tref, double smaller_range, double larger_range,
                                                 double period) {
    double prob = get_prob_smaller_time(op1, tref, smaller_range, period) * get_prob_larger_time(op2, tref, larger_range, period)
       + get_prob_larger_time(op1, tref, larger_range, period) * get_prob_smaller_time(op2, tref, smaller_range, period);
    return prob / 2.0;
  }

  // -----------------------------------

  // Draw a random time from the joint pdf defined through split_times the double_pdf function
  inline double get_close_time(triqs::mc_tools::random_generator &rng, fop_t const &op, std::vector<fop_t> const &op_list, double tmax) {
    double tau_ref = op_list[rng(op_list.size())].tau;
    return wrap(tau_ref + double_icdf(rng(), op.left_width, op.right_width, tmax), tmax);
  }

  // Draw a random time larger than tref using the pdf
  inline double get_close_larger_time(triqs::mc_tools::random_generator &rng, fop_t const &op, double tref, double range, double period) {
    return wrap(tref + icdf(rng(), op.left_width, range), period);
  }

  // Draw a random time smaller than tref using the pdf
  inline double get_close_smaller_time(triqs::mc_tools::random_generator &rng, fop_t const &op, double tref, double range, double period) {
    return wrap(tref - icdf(rng(), op.right_width, range), period);
  }

  // Draw a random time smaller and larger than a randomly chosen split_time
  inline std::array<double, 2> get_close_smaller_and_larger_time(triqs::mc_tools::random_generator &rng, fop_t const &op1, fop_t const &op2,
                                                                     std::vector<double> const &split_times, double tmax) {
    double tref = split_times[rng(split_times.size())];
    if (rng(2))
      return {get_close_smaller_time(rng, op1, tref, tmax / 2.0, tmax), get_close_larger_time(rng, op2, tref, tmax / 2.0, tmax)};
    else
      return {get_close_larger_time(rng, op1, tref, tmax / 2.0, tmax), get_close_smaller_time(rng, op2, tref, tmax / 2.0, tmax)};
  }

  // Draw a random time smaller and larger than a randomly chosen split_time
  inline std::array<double, 2> get_close_smaller_and_larger_time(triqs::mc_tools::random_generator &rng, fop_t const &op1, fop_t const &op2,
                                                                     double tref, double smaller_range, double larger_range, double period) {
    if (rng(2))
      return {get_close_smaller_time(rng, op1, tref, smaller_range, period), get_close_larger_time(rng, op2, tref, larger_range, period)};
    else
      return {get_close_larger_time(rng, op1, tref, larger_range, period), get_close_smaller_time(rng, op2, tref, smaller_range, period)};
  }

  // -----------------------------------

  // Draw a random time from the joint pdf defined through split_times the double_pdf function
  // and return both the time and its proposition probability
  inline std::array<double, 2> get_close_time_and_prob(triqs::mc_tools::random_generator &rng, fop_t op, std::vector<fop_t> const &op_list,
                                                           double tmax) {
    op.tau = get_close_time(rng, op, op_list, tmax);
    return {op.tau, get_prob(op, op_list, tmax)};
  }

  // Draw a random time larger than tref using the pdf
  // and return both the time and its proposition probability
  inline std::array<double, 2> get_close_larger_time_and_prob(triqs::mc_tools::random_generator &rng, fop_t op, double tref, double range,
                                                                  double period) {
    op.tau = get_close_larger_time(rng, op, tref, range, period);
    return {op.tau, get_prob_larger_time(op, tref, range, period)};
  }

  // Draw a random time smaller than tref using the pdf
  // and return both the time and its proposition probability
  inline std::array<double, 2> get_close_smaller_time_and_prob(triqs::mc_tools::random_generator &rng, fop_t op, double tref, double range,
                                                                   double period) {
    op.tau = get_close_smaller_time(rng, op, tref, range, period);
    return {op.tau, get_prob_smaller_time(op, tref, range, period)};
  }

  // Draw a random time smaller and larger than a randomly chosen split_time
  // and return both the times and the proposition probability
  inline std::array<double, 3> get_close_smaller_and_larger_time_and_prob(triqs::mc_tools::random_generator &rng, fop_t op1, fop_t op2,
                                                                                       std::vector<double> const &split_times, double tmax) {
    auto [t1, t2] = get_close_smaller_and_larger_time(rng, op1, op2, split_times, tmax);
    op1.tau       = t1;
    op2.tau       = t2;
    return {t1, t2, get_prob_smaller_and_larger_time(op1, op2, split_times, tmax)};
  }

  // Draw a random time smaller and larger than a randomly chosen split_time
  // and return both the times and the proposition probability
  inline std::array<double, 3> get_close_smaller_and_larger_time_and_prob(triqs::mc_tools::random_generator &rng, fop_t op1, fop_t op2,
                                                                                       double tref, double smaller_range, double larger_range,
                                                                                       double period) {
    auto [t1, t2] = get_close_smaller_and_larger_time(rng, op1, op2, tref, smaller_range, larger_range, period);
    op1.tau       = t1;
    op2.tau       = t2;
    return {t1, t2, get_prob_smaller_and_larger_time(op1, op2, tref, smaller_range, larger_range, period)};
  }

} // namespace inchworm
