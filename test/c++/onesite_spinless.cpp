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

TEST(inchworm, Hubbard_1site_spinless) {

  constr_params_t cp;
  cp.beta      = 2.0;
  cp.gf_struct = {{"up", {0}}};
  cp.n_tau     = 4;
  cp.n_iw      = 4;

  //mat_t theta   = {{1.0}};
  //vec_t epsilon = {3.0};
  mat_t theta = {{1.5, -1.0, 1.7}};
  vec_t epsilon = {-2.0, 0.4, 1.5};

  //n_site, n_bath, n_spin, U, mu, t
  auto [S, sp, u_tau] = test_setup(1, 3, 1, 0.0, 0.0, 0.0, cp, theta, epsilon);

  double tau_max   = cp.beta;
  double tau_split = cp.beta * 0.9;

  //solve_cthyb(S, sp, u_tau, tau_max);
  //solve_selfconsistent(S, sp, u_tau, tau_split, tau_max);
  solve_green(S, sp, u_tau);
}

MAKE_MAIN
