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
#include <triqs/atom_diag/partial_trace.hpp>

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

void print_eigensystems(triqs::atom_diag::atom_diag<false> const &ad) {
  for (auto sp : ad.get_eigensystems()) {
    for (auto l : sp.eigenvalues) { std::printf("% 2.3f ", l); }
    std::printf("\n\n");

    for (int i = 0; i < sp.eigenvalues.size(); i++) {
      for (int j = 0; j < sp.eigenvalues.size(); j++) { std::printf("% 2.3f ", sp.unitary_matrix(i, j)); }
      std::printf("\n");
    }
    //for (auto u : sp.unitary_matrix) { TRIQS_PRINT(u); }
    std::printf("\n\n");
  }
  std::printf("\n");
}

void print_matrix(triqs::arrays::matrix<double> m) {
  for (int i = 0; i < first_dim(m); i++) {
    for (int j = 0; j < second_dim(m); j++) { std::printf("% 2.3f ", m(i, j)); }
    std::printf("\n");
  }
  std::printf("\n\n");
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

TEST(inchworm, partial_trace) {

  //int n_bath             = 1;
  //double theta[n_bath]   = {0.5};
  //double epsilon[n_bath] = {0.1};
  int n_bath       = 1;
  double theta[]   = {0.5, 0.5, 0.5};
  double epsilon[] = {0.1, 0.1, 0.1};
  double mu        = 0.0;
  double U         = 1.0;
  auto fops        = make_fops(n_bath + 1);
  double dtau      = 1.0;

  auto h = U * (n("up", 0) * (n("dn", 0))); // 0 is the only interacting orbital
  h -= mu * (n("up", 0) + n("dn", 0));

  for (int i = 0; i < n_bath; i++) {
    h += 10000*theta[i] * (c_dag("up", 0) * c("up", i + 1) + c_dag("up", i + 1) * c("up", 0));
    h += 10000*theta[i] * (c_dag("dn", 0) * c("dn", i + 1) + c_dag("dn", i + 1) * c("dn", 0));
    h += epsilon[i] * (n("up", i + 1) + n("dn", i + 1));
    h -= mu * (n("up", i + 1) + n("dn", i + 1));
  }

  std::vector<many_body_op_t> qn;
  qn.resize(1);
  qn[0] += n("up", 0) + n("dn", 0);

  auto ad = triqs::atom_diag::atom_diag<false>(h, fops, qn);
  auto es = ad.get_eigensystems();

  //print_eigensystems(ad);
  /*
  triqs::arrays::matrix<double> expected_partial_sum(4, 4);
  expected_partial_sum = 0;
  for (int i = 0; i < 4; i++) expected_partial_sum(i, i) = epsilon[0];
  expected_partial_sum(3, 3) += U;
  print_matrix(expected_partial_sum);
  */

  //auto ps = partial_sum(ad, 2, [](double x) { return x; });
  //print_matrix(ps);
  //std::cout << "Are they equal?\n" << (expected_partial_sum != ps) << "\n"; // this is weird, this should be equal, something I do not get for now

  auto ps = partial_sum(ad, 2, [dtau](double x) { return std::exp(-dtau * x); });
  print_matrix(ps);
}

MAKE_MAIN
