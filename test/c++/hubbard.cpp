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
#include <inchworm/util.hpp>

#include <triqs/gfs.hpp>
#include <triqs/h5.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace inchworm;
//namespace nda = triqs::arrays;
using mat_t = triqs::arrays::array<double, 2>;
using vec_t = triqs::arrays::array<double, 1>;

// Prepare funcdamental operator set
std::pair<fundamental_operator_set, std::vector<many_body_op_t>> make_fops(int n_site, int n_bath, int linear_index, int n_spin) {
  fundamental_operator_set fops;
  std::vector<many_body_op_t> qn;
  qn.resize(1);
  for (int spin = 0; spin < n_spin; spin++)
    for (int i = 0; i < n_site; i++) {
      auto sp = ((spin == 0) ? "up" : "dn");
      fops.insert(sp, i);
      qn[0] += n(sp, i);
    }
  for (int spin = 0; spin < n_spin; spin++)
    for (int i = linear_index; i < linear_index + n_bath; i++) {
      auto sp = ((spin == 0) ? "up" : "dn");
      fops.insert(sp, i);
      qn[0] += n(sp, i);
    }
  return std::pair<fundamental_operator_set, std::vector<many_body_op_t>>(fops, qn);
}

void self_consistent_hubbard(int n_site, int n_bath, int n_spin, double U, double mu, double t, constr_params_t const &cp, mat_t const &theta,
                             vec_t const &epsilon, double tau_max, double tau_split) {
  // Set up the Solver
  solver_core S(cp);

  for (auto const &tau : S.Delta_tau[0].mesh()) {
    double val;

    for (int block = 0; block < cp.gf_struct.size(); block++) {
      S.Delta_tau[block][tau] = 0.0;
      for (int i = 0; i < n_site; i++) {
        for (int j = 0; j < n_site; j++) {
          for (int n = 0; n < n_bath; n++) {
            if (epsilon(n) >= 0.0) // to avoid numerical instability, assign hyb differently depending on the sign of epsilon(n).
              val = -theta(i, n) * theta(j, n) * (std::exp(-((double)tau) * (epsilon(n))) / (1. + std::exp(-cp.beta * epsilon(n))));
            else
              val = -theta(i, n) * theta(j, n) * (std::exp(-((double)tau - cp.beta) * (epsilon(n))) / (1. + std::exp(cp.beta * epsilon(n))));
            S.Delta_tau[block][tau](i, j) += val;
          }
        }
      }
    }
  }

  for (int block = 0; block < cp.gf_struct.size(); block++)
    for (int i = 0; i < n_site; i++)
      for (int j = 0; j < n_site; j++) std::printf("%d %d % 4.8f\n", i, j, S.Delta_tau[block][cp.n_tau - 1](i, j));
  //exit(0);

  auto h_atom = 0 * n("up", 0);
  for (int j = 0; j < n_site; j++) {
    h_atom -= mu * n("up", j);

    if (n_spin == 2) {
      h_atom -= mu * n("dn", j);
      h_atom += U * n("up", j) * n("dn", j);
    }
    for (int i = 0; i < n_site; i++) {
      if (i != j) {
        h_atom -= t * c_dag("up", i) * c("up", j);
        if (n_spin == 2) h_atom -= t * c_dag("dn", i) * c("dn", j);
      }
    }
  }

  // Compare against the reference data
  // h5diff("hubbard.out.h5", "hubbard.ref.h5")
  auto [fops_tot, qn_tot]   = make_fops(n_site, n_bath, n_site, n_spin);
  auto [fops_atom, qn_atom] = make_fops(n_site, 0, n_site, n_spin);
  auto [fops_bath, qn_bath] = make_fops(0, n_bath, n_site, n_spin);

  // Solve Parameters
  solve_params_t sp;
  sp.h_int           = h_atom;
  sp.n_cycles        = 50000;
  sp.length_cycle    = 10;
  sp.n_warmup_cycles = 20;
  sp.max_time        = -1;
  sp.verbosity       = 3;
  sp.post_process    = true;
  sp.measure_sign    = true;
  sp.quantum_numbers = qn_atom;
  sp.random_seed     = 12345789 + 928374 * mpi::communicator().rank();

  auto h_hyb  = 0.0 * n("up", 0);
  auto h_bath = 0.0 * n("up", n_site);

  for (int i = 0; i < n_site; i++) {
    for (int k = 0; k < n_bath; k++) {
      h_hyb += theta(i, k) * (c_dag("up", i) * c("up", k + n_site));
      h_hyb += theta(i, k) * (c_dag("up", k + n_site) * c("up", i));
      if (n_spin == 2) {
        h_hyb += theta(i, k) * (c_dag("dn", i) * c("dn", k + n_site));
        h_hyb += theta(i, k) * (c_dag("dn", k + n_site) * c("dn", i));
      }
    }
  }

  for (int k = 0; k < n_bath; k++) {
    h_bath += epsilon(k) * n("up", k + n_site);
    if (n_spin == 2) h_bath += epsilon(k) * n("dn", k + n_site);
  }

  auto ad_tot  = triqs::atom_diag::atom_diag<false>(h_atom + h_bath + h_hyb, fops_tot);
  auto ad_atom = triqs::atom_diag::atom_diag<false>(h_atom, fops_atom, qn_atom);
  auto ad_bath = triqs::atom_diag::atom_diag<false>(h_bath, fops_bath);

  u_tau_t u_tau = make_ED_propagator(ad_tot, ad_atom, ad_bath, cp.beta, cp.n_tau);
  std::printf("\n##################\nexact U(beta):\n");
  print(u_tau, tau_max);

  // Solve the impurity model
  auto result_cthyb = S.solve_cthyb(sp, tau_max);
  auto result_sc    = S.solve_self_consistently(sp, u_tau, tau_split, tau_max);

  for (int bl = 0; bl < result_sc.u_frame.size(); bl++) EXPECT_ARRAY_NEAR(((matrix_t)u_tau[bl][cp.n_tau - 1]), result_cthyb.u_frame[bl], 0.05*u_tau[0][cp.n_tau - 1](0,0)); // U[0](0,0) is essentially always the biggest value
  for (int bl = 0; bl < result_sc.u_frame.size(); bl++) EXPECT_ARRAY_NEAR(((matrix_t)u_tau[bl][cp.n_tau - 1]), result_sc.u_frame[bl], 0.05*u_tau[0][cp.n_tau - 1](0,0) );
}

//*
TEST(inchworm, Hubbard_1site_spinless) {

  constr_params_t cp;
  cp.beta      = 2.0;
  cp.gf_struct = {{"up", {0}}};
  cp.n_tau     = 500;
  cp.n_iw      = 250;
  //cp.n_step      = 250;

  mat_t theta = {{1.5, -1.0, 1.7}};
  //vec_t epsilon = {0.0, 0.0, 0.0};
  vec_t epsilon = {-2.0, 0.4, 1.5};
  self_consistent_hubbard(1, 3, 1, 0.0, 0.0, 0.0, cp, theta, epsilon, cp.beta, cp.beta * 0.9);
  //self_consistent_hubbard(1, 3, 1, 0.0, 2.0, 0.0, cp, theta, epsilon, cp.beta, cp.beta * 0.9);
}
//*/

//*
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
//*/

//*
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
//*/

// long:
//*
TEST(inchworm, Hubbard_2sites) { // NOLINT

  constr_params_t cp;
  cp.beta = 10.0;
  //cp.gf_struct = {{"up", {0, 1}}};
  cp.gf_struct = {{"up", {0, 1}}, {"dn", {0, 1}}};
  cp.n_tau     = 500;
  cp.n_iw      = 250;

  triqs::arrays::array<double, 2> theta   = {{0.1, -0.3, -0.4}, {0.1, 0.2, 0.4}};
  triqs::arrays::array<double, 1> epsilon = {1.0, -1.0, 1.2};
  //self_consistent_hubbard(2, 2, 2, 0.0, 0.0, 0.0, cp, theta, epsilon, cp.beta, cp.beta * 0.9);
  self_consistent_hubbard(2, 3, 2, 4.0, -3.0, 1.0, cp, theta, epsilon, cp.beta, cp.beta * 0.9);
}
//*/

MAKE_MAIN
