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
#include <inchworm/solver_core.hpp>

#include <triqs/gfs.hpp>
#include <triqs/h5.hpp>
#include <triqs/test_tools/gfs.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/arrays/blas_lapack/dot.hpp>

using namespace inchworm;

void print_energies(std::vector<std::vector<double>> const &E) {
  for (auto sp : E) {
    for (auto l : sp) { std::printf("% 2.3f ", l); }
    std::printf("\n");
  }
  std::printf("\n");
}

// Prepare funcdamental operator set
fundamental_operator_set make_fops(int N) {
  fundamental_operator_set fops;
  for (int o : range(N)) {
    fops.insert("up", o);
    fops.insert("dn", o);
  }
  return fops;
}

TEST(atom_diag_real, atom_diag) {

  int n_bath             = 2;
  double theta[n_bath]   = {0.5, 0.8};
  double epsilon[n_bath] = {0.1, -0.1};
  double mu              = 0.0;
  double U               = 8.0;
  auto fops              = make_fops(n_bath + 1);

  auto h = U * (n("up", 0) * (n("dn", 0))); // 0 is the only interacting orbital

  for (int i = 0; i < n_bath; i++) {
    h += theta[i] * (c_dag("up", 0) * c("up", i + 1) + c_dag("up", i + 1) * c("up", 0));
    h += theta[i] * (c_dag("dn", 0) * c("dn", i + 1) + c_dag("dn", i + 1) * c("dn", 0));
    h += epsilon[i] * (n("up", i + 1) + n("dn", i + 1));
  }

  PRINT_TRIQS(h);
  /*
  auto ad = triqs::atom_diag::atom_diag<false>(h, fops);

  auto vac = ad.get_vacuum_state();
  EXPECT_NEAR(1, dot(vac, vac), 1e-14);

  std::vector<matrix<double>> V      = ad.get_unitary_matrices();
  std::vector<std::vector<double>> E = ad.get_energies();

  print_energies(E);
  */
}

MAKE_MAIN;
