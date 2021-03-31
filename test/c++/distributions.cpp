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

#include <inchworm/distributions.hpp>

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
        for (auto w2 : wlst)
	  EXPECT_NEAR(double_pdf(x * tmax, w1, w2, tmax), 1.0 / tmax, tol(tmax, std::min(w1,w2)));
      }
    }
  }
}

TEST(Distributions, norm_pdf) {

  std::vector<double> tmaxlst = {1.0, 10.0, 100.0};
  std::vector<double> wlst    = {0.1, 10.0, 1000.0};

  for (auto tmax : tmaxlst) {
    for (auto w1 : wlst) {

      long Nint = 1e+6;
      double dtau = tmax / Nint;

      double int_pdf        = 0.0;
      for (auto i: range(Nint)) int_pdf += pdf((0.5 + i) * dtau, w1, tmax) * dtau;
      EXPECT_NEAR(int_pdf, 1.0, 1.0/Nint);

      for (auto w2 : wlst) {
        double int_double_pdf = 0.0;
        for (auto i: range(Nint)) int_double_pdf += double_pdf((0.5 + i) * dtau, w1, w2, tmax) * dtau;
        EXPECT_NEAR(int_double_pdf, 1.0, 1.0/Nint);
      }
    }
  }
}

TEST(Distributions, icdf) {

  std::vector<double> tmaxlst = {1.0, 10.0};
  std::vector<double> wlst    = {0.1, 10.0, 100.0};

  for (auto tmax : tmaxlst) {
    for (auto w1 : wlst) {

      double tol = 1e-14;

      EXPECT_NEAR(icdf(0.0, w1, tmax), 0.0, tol);
      EXPECT_NEAR(icdf(1.0, w1, tmax), tmax, tol);

      for (auto w2 : wlst) {
	EXPECT_NEAR(double_icdf(0.0, w1, w2, tmax), 0.0, tol);
	EXPECT_NEAR(double_icdf(0.5, w1, w2, tmax), tmax/2.0, tol);
	EXPECT_NEAR(double_icdf(1.0, w1, w2, tmax), tmax, tol);
      }
    }
  }
}
