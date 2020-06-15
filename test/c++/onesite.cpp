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

TEST(inchworm, Hubbard_1site) {

  constr_params_t cp;
  cp.beta      = 2.0;
  cp.gf_struct = {{"up", {0}}, {"dn", {0}}};
  cp.n_tau     = 500;
  cp.n_iw      = 250;

  mat_t theta   = {{0.9, -1.0, 1.1}};
  vec_t epsilon = {1.0, -2.0, 0.0};
  self_consistent_hubbard(1, 3, 2, 4.0, -2.0, 0.0, cp, theta, epsilon, cp.beta, cp.beta * 0.9);
}

MAKE_MAIN
