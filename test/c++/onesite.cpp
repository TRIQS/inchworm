/*******************************************************************************
 *
 * inchworm: A TRIQS based impurity solver
 *
 * Copyright (c) 2019 The Simons foundation
 *   authors: Maxime Charlebois
 *
 * inchworm is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * inchworm is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * inchworm. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "./hubbard.hpp"

TEST(inchworm, Hubbard_1site) { // NOLINT

  constr_params_t cp;
  cp.beta        = 1.0;
  cp.gf_struct   = {{"up", 1}, {"dn", 1}};
  cp.n_tau_green = 5;
  cp.n_tau_inch  = 21;
  cp.n_tau       = 10001;

  mat_t theta   = {{0.8, 1.0, 1.4}};
  vec_t epsilon = {-2.0, 0.0, 4.0};

  int n_site = 1;
  int n_bath = epsilon.size();
  int n_spin = cp.gf_struct.size();
  double U   = 1.0;
  double mu  = 1.0;
  double t   = 1.0;

  auto [S, sp, u_tau, G_tau] = test_setup(n_site, n_bath, n_spin, U, mu, t, cp, theta, epsilon);

  double tau_max   = cp.beta;
  double tau_split = cp.beta * 0.9;

  // Test cthyb
  auto result_cthyb = S.solve_cthyb(sp, tau_max);
  EXPECT_TRUE(relative_distance(get_frame(u_tau, cp.n_tau_inch -1), result_cthyb.frame) < 0.05);

  // Test selfconsistent
  auto result_sc = S.solve_self_consistently(sp, u_tau, tau_split, tau_max);
  EXPECT_TRUE(relative_distance(get_frame(u_tau, cp.n_tau_inch -1), result_sc.frame) < 0.05);

  // Test inchworm
  S.solve_inchworm(sp);
  EXPECT_TRUE(relative_distance(S.u_tau, u_tau) < 0.05);

  // Test green
  S.u_tau = u_tau;
  S.solve_green(sp);
  EXPECT_BLOCK_GF_NEAR(S.G_tau, G_tau, 0.01);
}

MAKE_MAIN
