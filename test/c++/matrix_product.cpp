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

//#include <inchworm/solver_core.hpp>
#include <inchworm/mc/impurity_product.hpp>
#include <inchworm/types.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace inchworm;

void print_bin(int v) { std::cout << std::bitset<4>(v); }

void print_ad(triqs::atom_diag::atom_diag<false> const &ad) {
  printf("\n%d \n", ad.get_full_hilbert_space_dim());
  for (auto fock_sp : ad.get_fock_states()) {
    for (auto l : fock_sp) {
      print_bin(l);
      std::printf(" ");
    }
    std::printf("\n");
  }
  std::printf("\n\n");
}

TEST(inchworm, matrix_product) {

  int n_times = 4;
  //  double beta = 1.0;
  double mu = 0.0;
  double U  = 8.0;
  double t  = 1.0;
  fundamental_operator_set fops;
  int n_site = 2;

  auto qn_vector = std::vector<triqs::operators::many_body_operator_generic<scalar_t>>();
  auto h         = 0 * (n("up", 0));
  auto n_tot     = 0 * (n("up", 0));

  for (int i = 0; i < n_site; i++) {
    fops.insert("up", i);
    fops.insert("dn", i);
    h += U * (n("up", i) * (n("dn", i)));
    h -= mu * (n("up", i) + n("dn", i));
    for (int j = 0; j < n_site; j++) {
      if (i != j) h -= t * (c_dag("up", i) * (c("up", j)) + c_dag("dn", i) * (c("dn", j)));
    }
    n_tot += n("up", i) + n("dn", i);
  }
  qn_vector.push_back(n_tot);

  inchworm::atom_diag ad = {h, fops, qn_vector};
  //auto propagator_struct = find_propagator_struct(ad);
  //auto u_tau             = u_tau_t{{beta, Fermion, n_times}, propagator_struct};
  u_tau_t u_tau = make_propagator(ad, n_times);
  //for (int i = 0; i < ad.n_subspaces(); i++) { std::cout << propagator[i] << "\n"; }
  //print_block_gf_first_time(propagator);

  //print_ad(ad);
  //print_block_gf_first_time(propagator);
  std::vector<diagram::time_and_index_t> c    = {{0.00, 0}, {0.05, 0}};
  std::vector<diagram::time_and_index_t> cdag = {{0.01, 0}, {0.08, 0}};
  std::vector<double> split_times             = {};
  time_diagram_t diagram(c, cdag, split_times);

  u_frame_t u_frame = make_zero_propagator_frame(ad);
  u_frame           = propagator_product(ad, diagram, 0.1);
  for (int bl = 0; bl < u_frame.size(); bl++) u_tau[bl][1] = u_frame[bl];
  std::cout << u_frame;
  std::cout << frobenius_norm(u_frame) << "\n";

  u_frame = propagator_product(ad, diagram, 0.2, &u_tau);
  for (int bl = 0; bl < u_frame.size(); bl++) u_tau[bl][2] = u_frame[bl];
  //assign_frame_to_propagator(u_tau, u_frame, 2);
  std::cout << u_frame;
  std::cout << frobenius_norm(u_frame) << "\n";

  u_frame = propagator_product(ad, diagram, 0.3, &u_tau);
  for (int bl = 0; bl < u_frame.size(); bl++) u_tau[bl][3] = u_frame[bl];
  //assign_frame_to_propagator(u_tau, u_frame, 3);
  std::cout << u_frame;
  std::cout << frobenius_norm(u_frame) << "\n";

  //propagator_product(8.0, ad, diagram, 0.4); //////////ATTENTION CA MARCHE (POURQUOI?)
}

MAKE_MAIN
