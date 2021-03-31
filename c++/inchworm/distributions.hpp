#pragma once

#include "./types.hpp"

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

} // namespace inchworm
