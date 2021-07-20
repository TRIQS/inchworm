/*******************************************************************************
 *
 * inchworm: A TRIQS based impurity solver
 *
 * Copyright (c) 2019 The Simons foundation
 *   authors: Nils Wentzell
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

#include <inchworm/types.hpp>
#include <inchworm/u_frame.hpp>
#include <inchworm/interpolator.hpp>

#include <gtest/gtest.h>

using namespace inchworm;

void test_function(auto &&f, double tol = 0.001) {

  auto gf_struct = gf_struct_t{{"up", 2}, {"dn", 2}};

  // Setup BlockGf
  double beta = 10.0;
  int n_tau   = 101;
  auto u_tau  = u_tau_t{{beta, Fermion, n_tau}, gf_struct};

  // Init BlockGf & Interpolator
  u_tau[bl_][tau_](i_, j_) << f(bl_, tau_, i_, j_);
  auto u_interpolator = interpolator_t<scalar_t>{u_tau, n_tau};

  auto tau_mesh_fine = triqs::mesh::imtime{beta, Fermion, 10 * n_tau};
  auto u_tau_interp  = u_tau_t{tau_mesh_fine, gf_struct};
  auto u_tau_exact   = u_tau_t{tau_mesh_fine, gf_struct};

  // Initialize Green Functions on Finer Mesh
  // using both interpolation and the exact function
  for (auto const &tau : tau_mesh_fine)
    for (auto bl : range(gf_struct.size())) {
      auto n_orb = gf_struct[bl].second;
      for (auto [i, j] : product_range(n_orb, n_orb)) {
        u_tau_interp[bl][tau](i, j) = u_interpolator(bl, tau, i, j);
        u_tau_exact[bl][tau](i, j)  = f(bl, tau, i, j);
      }
    }
  EXPECT_LT(relative_distance(u_tau_interp, u_tau_exact), tol);

  //{
  //auto f = h5::file{"interp.h5", 'w'};
  //h5_write(f, "u_interp", u_tau_interp);
  //h5_write(f, "u_exact", u_tau_exact);
  //}
}

TEST(Interpolator, Lin) {
  // Linear Function
  auto f_lin = [](auto bl, auto tau, auto i, auto j) { return tau * (bl + i + j); };
  test_function(f_lin, 1e-14);
}

TEST(Interpolator, Quad) {
  // Quadratic Function
  auto f_quad = [](auto bl, auto tau, auto i, auto j) { return (tau + 1.0) * (tau + 1.0) * (bl + i + j); };
  test_function(f_quad, 0.002);
}

TEST(Interpolator, Cubic) {
  // Cubic function
  auto f_cubic = [](auto bl, auto tau, auto i, auto j) { return (tau + 1.0) * (tau + 1.0) * (tau + 1.0) * (bl + i + j); };
  test_function(f_cubic, 0.01);
}

TEST(Interpolator, Exp) {
  // Exponential function
  auto f_exp = [](auto bl, auto tau, auto i, auto j) { return exp(tau * (bl + i + j)); };
  test_function(f_exp, 0.01);
}
