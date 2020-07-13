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
  cp.n_tau     = 4;
  cp.n_iw      = 5;

  //mat_t theta   = {{0.9, -1.0, 1.1}};
  //vec_t epsilon = {1.0, -2.0, 0.0};

  mat_t theta   = {{1.0}};
  vec_t epsilon = {1.0};

  //n_site, n_bath, n_spin, U, mu, t
  int n_site = 1; 
  int n_bath = epsilon.size();
  int n_spin = cp.gf_struct.size();
  double U = 0.0;
  double mu = 1.0; 
  double t = 1.0;

  auto [S, sp, u_tau] = test_setup(n_site, n_bath, n_spin, U, mu, t, cp, theta, epsilon);

  double tau_max   = cp.beta;
  double tau_split = cp.beta * 0.9;

  //solve_cthyb(S, sp, u_tau, tau_max);
  //solve_selfconsistent(S, sp, u_tau, tau_split, tau_max);
  S.solve_green(sp, u_tau);

  auto G_tau_exact = green_U0_setup(n_site, n_bath, n_spin, mu, t, cp, theta, epsilon);

  for(auto const & tau: G_tau_exact[0].mesh()){
    std::cout << "\n\nCalculated: " << S.G_tau[0][tau] << "\nExact: " << G_tau_exact[0][tau];
  }

}

MAKE_MAIN
