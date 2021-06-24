/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2020 Simons Foundation
 *   author: N. Wentzell
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#undef NDEBUG

#include <nda/gtest_tools.hpp>

#include <triqs/stat/histograms.hpp>

#include <inchworm/distributions.hpp>

#include <random>

using namespace inchworm;

TEST(Distributions, pdf) {

  std::vector<double> tmaxlst = {1.0, 10.0, 100.0};
  std::vector<double> xlst    = {0.0, 0.3, 0.7, 1.0};
  std::vector<double> wlst    = {1e+2, 1e+4, 1e+8};

  auto tol = [](double tmax, double w) { return 1 - std::exp(-tmax / w); };

  for (auto x : xlst) {
    for (auto tmax : tmaxlst) {
      for (auto w1 : wlst) {
        EXPECT_NEAR(pdf(x * tmax, w1, tmax), 1.0 / tmax, tol(tmax, w1));
        for (auto w2 : wlst) EXPECT_NEAR(double_pdf(x * tmax, w1, w2, tmax), 1.0 / tmax, tol(tmax, std::min(w1, w2)));
      }
    }
  }
}

double integrate(auto f, double tmin, double tmax, long Nint = 1e+6) {
  double dtau    = (tmax - tmin) / Nint;
  double int_pdf = 0.0;
  for (auto i : range(Nint)) int_pdf += f(tmin + (0.5 + i) * dtau) * dtau;
  return int_pdf;
}

TEST(Distributions, norm_pdf) {

  std::vector<double> tmaxlst = {1.0, 10.0, 100.0};
  std::vector<double> wlst    = {0.1, 10.0, 1000.0, 1e+8};

  for (auto tmax : tmaxlst) {
    for (auto w1 : wlst) {

      auto f         = [&](double tau) { return pdf(tau, w1, tmax); };
      double int_pdf = integrate(f, 0, tmax);
      EXPECT_NEAR(int_pdf, 1.0, 1.0 / 1e+6);

      for (auto w2 : wlst) {
        auto g                = [&](double tau) { return double_pdf(tau, w1, w2, tmax); };
        double int_double_pdf = integrate(g, 0, tmax);
        EXPECT_NEAR(int_double_pdf, 1.0, 1.0 / 1e+6);
      }
    }
  }
}

TEST(Distributions, icdf) {

  std::vector<double> tmaxlst = {1.0, 10.0, 100.0};
  std::vector<double> wlst    = {0.1, 10.0, 1000.0, 1e+6};

  for (auto tmax : tmaxlst) {
    for (auto w1 : wlst) {

      double tol = 1e-10;

      EXPECT_NEAR(icdf(0.0, w1, tmax), 0.0, tol);
      EXPECT_NEAR(icdf(1.0, w1, tmax), tmax, tol);

      for (auto w2 : wlst) {
        EXPECT_NEAR(double_icdf(0.0, w1, w2, tmax), 0.0, tol);
        EXPECT_NEAR(double_icdf(0.5, w1, w2, tmax), tmax / 2.0, tol);
        EXPECT_NEAR(double_icdf(1.0, w1, w2, tmax), tmax, tol);
      }
    }
  }
}

TEST(Distributions, icdf_inverse) {

  std::vector<double> tmaxlst = {1.0, 10.0, 100.0};
  std::vector<double> wlst    = {0.1, 10.0, 1000.0, 1e+6};

  for (auto tmax : tmaxlst) {
    for (auto w1 : wlst) {

      // -- Calculate the cumulative distrubition function numerically

      auto f   = [&](double tau) { return pdf(tau, w1, tmax); };
      auto cdf = [&](double tau) { return integrate(f, 0, tau); };

      // -- Check that inverse holds both directions

      double tol = 1e-6;
      EXPECT_NEAR((cdf(icdf(0.2, w1, tmax))), 0.2, tol);
      EXPECT_NEAR((cdf(icdf(0.5, w1, tmax))), 0.5, tol);
      EXPECT_NEAR((cdf(icdf(0.7, w1, tmax))), 0.7, tol);

      tol = 1e-6 * tmax;
      // CAUTION: Do not evaluate deep in the exponential tail of the pdf
      if (w1 < tmax)
        EXPECT_NEAR((icdf(cdf(w1), w1, tmax)), w1, tol);
      else
        EXPECT_NEAR((icdf(cdf(tmax), w1, tmax)), tmax, tol);

      if (5.0 * w1 < tmax)
        EXPECT_NEAR((icdf(cdf(5.0 * w1), w1, tmax)), 5.0 * w1, tol);
      else
        EXPECT_NEAR((icdf(cdf(0.2 * tmax), w1, tmax)), 0.2 * tmax, tol);

      if (10.0 * w1 < tmax)
        EXPECT_NEAR((icdf(cdf(10.0 * w1), w1, tmax)), 10.0 * w1, tol);
      else
        EXPECT_NEAR((icdf(cdf(0.1 * tmax), w1, tmax)), 0.1 * tmax, tol);

      for (auto w2 : wlst) {

        // -- Calculate the cumulative distrubition function numerically

        auto g          = [&](double tau) { return double_pdf(tau, w1, w2, tmax); };
        auto double_cdf = [&](double tau) { return integrate(g, 0, tau); };

        // -- Check that inverse holds both directions

        tol = 1e-6;
        EXPECT_NEAR((double_cdf(double_icdf(0.2, w1, w2, tmax))), 0.2, tol);
        EXPECT_NEAR((double_cdf(double_icdf(0.5, w1, w2, tmax))), 0.5, tol);
        EXPECT_NEAR((double_cdf(double_icdf(0.7, w1, w2, tmax))), 0.7, tol);

        tol = 1e-6 * tmax;
        if (w1 < tmax / 2.0)
          EXPECT_NEAR((double_icdf(double_cdf(w1), w1, w2, tmax)), w1, tol);
        else
          EXPECT_NEAR((double_icdf(double_cdf(0.4 * tmax), w1, w2, tmax)), 0.4 * tmax, tol);

        if (w2 < tmax / 2.0)
          EXPECT_NEAR((double_icdf(double_cdf(tmax - w2), w1, w2, tmax)), tmax - w2, tol);
        else
          EXPECT_NEAR((double_icdf(double_cdf(0.6 * tmax), w1, w2, tmax)), 0.6 * tmax, tol);

        if (5.0 * w1 < tmax / 2.0)
          EXPECT_NEAR((double_icdf(double_cdf(5.0 * w1), w1, w2, tmax)), 5.0 * w1, tol);
        else
          EXPECT_NEAR((double_icdf(double_cdf(0.1 * tmax), w1, w2, tmax)), 0.1 * tmax, tol);

        if (5.0 * w2 < tmax / 2.0)
          EXPECT_NEAR((double_icdf(double_cdf(tmax - 5.0 * w2), w1, w2, tmax)), tmax - 5.0 * w2, tol);
        else
          EXPECT_NEAR((double_icdf(double_cdf(0.9 * tmax), w1, w2, tmax)), 0.9 * tmax, tol);
      }
    }
  }
}

TEST(Distributions, get_prob_norm) {

  std::vector<double> tmaxlst = {1.0, 10.0, 100.0};
  std::vector<double> wlst    = {0.1, 10.0, 1000.0, 1e+8};

  for (auto tmax : tmaxlst) {
    for (auto w1 : wlst) {

      auto d     = fop_t{0.0, false, 0, 0, 0, w1, w1};

      auto tau_splits = std::vector{0.4 * tmax, 0.7 * tmax, 0.9 * tmax};

      double tol = 1e-6;
      auto f     = [&](double tau) {
        d.tau = tau;
        return get_prob(d, tau_splits, tmax);
      };
      EXPECT_NEAR(integrate(f, 0.0, tmax), 1.0, tol);

      // Test 2d integration for get_prob_smaller_and_larger_time
      // Takes several minutes ..
      //auto i = [&](double tau) {
      //auto j = [&](double tau_dag) {
      //d.tau = tau;
      //d_dag.tau = tau_dag;
      //return get_prob_smaller_and_larger_time(d, d_dag, tau_splits, tmax);
      //};
      //return integrate(j, 0, tmax, 1e+4);
      //};
      //EXPECT_NEAR(integrate(i, 0, tmax, 1e+4), 1.0, 1e-3);
    }
  }
}

struct stdrng_t {
  //Standard mersenne_twister_engine seeded with std::random_device
  std::mt19937 gen{std::random_device{}()};

  double operator()(double b = 1.0) { return std::uniform_real_distribution<>(0, b)(gen); }
  double operator()(double a, double b) { return std::uniform_real_distribution<>(a, b)(gen); }
  std::integral auto operator()(std::integral auto i) { return std::uniform_int_distribution<>(0, i - 1)(gen); }
};

TEST(Distributions, get_close_time) {

  double tmax = 1.0;

  std::vector<double> split_times{0.0, 0.1 * tmax, 0.4 * tmax, 0.8 * tmax};

  double w1     = 0.1 * tmax * 500;
  double w2     = 0.2 * tmax * 500;
  auto make_fop = [&](double tau) { return fop_t{tau, false, 0, 0, 0, w1, w2}; };
  auto d        = make_fop(0.0);

  long Nsamples = 1e+7;
  long Nbins    = 101;

  auto rng = triqs::mc_tools::random_generator{}; // stdrng_t{};

  auto hist = triqs::stat::histogram(0.0, tmax, Nbins);
  for (auto i : range(Nsamples)) hist << get_close_time(rng, d, split_times, tmax); // double_icdf(rng(), w1, w2, tmax);

  auto hist_norm = pdf(hist);
  std::vector<double> dat_pdf(Nbins);
  auto hist_norm_dat = hist_norm.data();
  for (auto n : range(Nbins)) {
    auto bin_width = tmax / (Nbins - 1) * ((n == 0 or n == Nbins - 1) ? 0.5 : 1);
    double tau     = hist_norm.mesh_point(n);
    double prob;
    if (n == 0)
      prob = get_prob(make_fop(tau), split_times, tmax); // double_pdf(tau, w1, w2, tmax);
    else if (n == Nbins - 1)                             // Evaluate only for 'tmax - eps'
      prob = get_prob(make_fop(tau - 1e-10), split_times, tmax);
    else // Be sure to symmetrize around bin center
      prob = 0.5 * (get_prob(make_fop(tau - 1e-10), split_times, tmax) + get_prob(make_fop(tau + 1e-10), split_times, tmax));
    dat_pdf[n] = prob;
    EXPECT_NEAR(hist_norm.data()[n] / bin_width, prob, 5 * std::sqrt(double(Nbins) / Nsamples));
    hist_norm_dat[n] *= 1.0 / bin_width;
  }

  //{
  //auto f = h5::file{"dat.h5", 'w'};
  //h5::write(f, "pdf", dat_pdf);
  //h5::write(f, "hist", hist_norm_dat);
  //}
}
