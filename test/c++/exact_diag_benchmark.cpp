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

uint64_t get_MSB(uint64_t a, int shift) { return (a >> shift); }

uint64_t get_LSB(uint64_t a, int shift) { return (a % (1 << shift)); }

// index: all indices above this index will be traced out.
void print_partial_sum_propagator(triqs::atom_diag::atom_diag<false> const &ad, int index, double tau) {
  //printf("%d ", ad.get_full_hilbert_space_dim());
  //printf("%d \n", (1 << index));
  //exit(0);
  int dim_partial = (1 << index);
  int dim_full    = ad.get_full_hilbert_space_dim();
  int factor      = dim_full / dim_partial;
  EXPECTS(dim_partial < dim_full);
  EXPECTS(dim_full % dim_partial == 0);

  triqs::arrays::matrix<double> partial_sum(dim_partial, dim_partial);
  partial_sum = 0;

  auto es = ad.get_eigensystems();
  auto fs = ad.get_fock_states();
  EXPECTS(es.size() == fs.size());

  for (int s = 0; s < ad.n_subspaces(); s++) {
    EXPECTS(es[s].eigenvalues.size() == fs[s].size());
    int size   = ad.get_subspace_dim(s);
    auto EUdag = dagger(es[s].unitary_matrix);
    for (int i = 0; i < size; i++)
      for (int j = 0; j < size; j++) EUdag(i, j) *= std::exp(-tau * (es[s].eigenvalues[i] + ad.get_gs_energy()));
    auto H = es[s].unitary_matrix * EUdag;

    for (int i = 0; i < size; i++) {
      uint64_t traced_idx1    = get_MSB(fs[s][i], index);
      uint64_t preserved_idx1 = get_LSB(fs[s][i], index);

      for (int j = 0; j < size; j++) {
        uint64_t traced_idx2    = get_MSB(fs[s][j], index);
        uint64_t preserved_idx2 = get_LSB(fs[s][j], index);
        if (traced_idx1 == traced_idx2) { partial_sum(preserved_idx1, preserved_idx2) += H(i, j) / factor; }
      }
    }
  }

  for (int i = 0; i < dim_partial; i++) {
    for (int j = 0; j < dim_partial; j++) { std::printf("% 2.3f ", partial_sum(i, j)); }
    std::printf("\n");
  }
  std::printf("\n\n");

  //for (auto fock_sp : ad.get_fock_states()) {
  //  for (auto l : fock_sp) { std::printf("%2d ", l); }
  //  std::printf("\n");
  //}
  //std::printf("\n\n");
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

  //int n_bath             = 1;
  //double theta[n_bath]   = {0.5};
  //double epsilon[n_bath] = {0.1};
  int n_bath             = 3;
  double theta[n_bath]   = {0.5, 0.5, 0.5};
  double epsilon[n_bath] = {0.1, 0.1, 0.1};
  double mu = 0.0;
  double U  = 8.0;
  auto fops = make_fops(n_bath + 1);

  auto h = U * (n("up", 0) * (n("dn", 0))); // 0 is the only interacting orbital
  h -= mu * (n("up", 0) + n("dn", 0));

  for (int i = 0; i < n_bath; i++) {
    h += theta[i] * (c_dag("up", 0) * c("up", i + 1) + c_dag("up", i + 1) * c("up", 0));
    h += theta[i] * (c_dag("dn", 0) * c("dn", i + 1) + c_dag("dn", i + 1) * c("dn", 0));
    h += epsilon[i] * (n("up", i + 1) + n("dn", i + 1));
    h -= mu * (n("up", i + 1) + n("dn", i + 1));
  }

  auto ad = triqs::atom_diag::atom_diag<false>(h, fops);
  auto es = ad.get_eigensystems();
  print_partial_sum_propagator(ad, 2, 0.4);
}

MAKE_MAIN;
