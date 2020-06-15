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

TEST(inchworm, Hubbard_2sites_spinless) { // NOLINT

  constr_params_t cp;
  cp.beta      = 2.0;
  cp.gf_struct = {{"up", {0, 1}}};
  //cp.gf_struct = {{"up", {0, 1}}, {"dn", {0, 1}}};
  cp.n_tau = 500;
  cp.n_iw  = 250;

  triqs::arrays::array<double, 2> theta   = {{0.9, 0.5}, {0.3, 1.1}};
  triqs::arrays::array<double, 1> epsilon = {0.9, -0.3};
  self_consistent_hubbard(2, 2, 1, 4.0, -3.0, 1.0, cp, theta, epsilon, cp.beta, cp.beta * 0.9);
}

MAKE_MAIN
