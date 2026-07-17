// Copyright (c) 2020--present, The Simons Foundation
// This file is part of inchworm and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.


#include "./hubbard.hpp"

TEST(inchworm, Hubbard_2sites_spinless) { // NOLINT

  constr_params_t cp;
  cp.beta        = 2.0;
  cp.gf_struct   = {{"up", 2}};
  cp.n_tau_green = 5;
  cp.n_tau_inch  = 21;
  cp.n_tau       = 10001;

  mat_t theta   = {{0.9, 0.5}, {0.3, 1.1}};
  vec_t epsilon = {0.9, -0.3};

  int n_site = 2;
  int n_bath = epsilon.size();
  int n_spin = cp.gf_struct.size();
  double U   = 0.0;
  double mu  = 1.0;
  double t   = 1.0;

  auto [S, sp, u_tau, G_tau] = test_setup(n_site, n_bath, n_spin, U, mu, t, cp, theta, epsilon);

  double tau_max   = cp.beta;
  double tau_split = cp.beta * 0.9;

  // Test cthyb
  auto result_cthyb = S.solve_cthyb(sp, tau_max);
  EXPECT_LT(relative_distance(get_frame(u_tau, cp.n_tau_inch - 1), result_cthyb.frame), 0.03);

  // Test selfconsistent
  auto result_sc = S.solve_self_consistently(sp, u_tau, tau_split, tau_max);
  EXPECT_LT(relative_distance(get_frame(u_tau, cp.n_tau_inch - 1), result_sc.frame), 0.02);

  // Test inchworm
  S.solve_inchworm(sp);
  EXPECT_LT(relative_distance(S.u_tau, u_tau), 0.06);

  // Test green
  S.u_tau = u_tau;
  S.solve_green(sp);
  EXPECT_BLOCK_GF_NEAR(S.G_tau, G_tau, 0.01);
}

MAKE_MAIN
