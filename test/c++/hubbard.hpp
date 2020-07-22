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

#include <inchworm/solver_core.hpp>
#include <inchworm/util.hpp>

#include <triqs/gfs.hpp>
#include <triqs/h5.hpp>
#include <triqs/test_tools/gfs.hpp>
#include <triqs/arrays/linalg/eigenelements.hpp>

using namespace inchworm;
//namespace nda = triqs::arrays;
using mat_t = triqs::arrays::matrix<double>;
using vec_t = triqs::arrays::array<double, 1>;

// Prepare funcdamental operator set
inline std::pair<fundamental_operator_set, std::vector<many_body_op_t>> make_fops(int n_site, int n_bath, int bath_offset, int n_spin) {
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
    for (int i = bath_offset; i < bath_offset + n_bath; i++) {
      auto sp = ((spin == 0) ? "up" : "dn");
      fops.insert(sp, i);
      qn[0] += n(sp, i);
    }
  return std::pair<fundamental_operator_set, std::vector<many_body_op_t>>(fops, qn);
}

double one_fermion(double tau, double eps, double beta){
  if (eps >= 0){
    return -std::exp(-tau * eps) / (1. + std::exp(-beta * eps));
  } else {
    return -std::exp((beta - tau) * eps) / (1. + std::exp(beta * eps));
  }
}

inline g_tau_t green_U0_setup(int n_site, int n_bath, int n_spin, double mu, double t, constr_params_t const &cp, mat_t const &theta, vec_t const &epsilon) {

  auto G_tau = g_tau_t{{cp.beta, Fermion, cp.n_tau_green}, cp.gf_struct};

  int n = n_site + n_bath;
  auto H = mat_t(n, n);
  H = 0.;

  // Impurity Energies
  for(int i = 0; i < n_site; ++i)
    H(i, i) += -mu;

  // Bath Energies
  for(int i = 0; i < n_bath; ++i)
    H(i + n_site, i + n_site) += epsilon(i);

  // Impurity Hopping
  for(int i = 0; i < n_site; ++i)
    for(int j = 0; j < n_site; ++j)
      if(i != j)
        H(i, j) += -t;

  // Coupling
  for(int i = 0; i < n_site; ++i)
    for(int j = 0; j < n_bath; ++j)
    {
        H(i, n_site + j) += theta(i, j);
        H(n_site + j, i) += theta(i, j);
    }

  auto [evals, evecs] = triqs::arrays::linalg::eigenelements(H);

  auto get_G_tau = [evals = evals, evecs = evecs, beta=cp.beta, n, n_site](double tau){
    auto G_full_diag = mat_t(n, n);
    G_full_diag = 0.;
    for(int i = 0; i < n; ++i)
      G_full_diag(i, i) = one_fermion(tau, evals(i), beta);
    auto G_full = mat_t{dagger(evecs) * G_full_diag * evecs};
    return mat_t{G_full(range(n_site), range(n_site))};
  };

  for(auto const & tau: G_tau[0].mesh()){
    for(int sp = 0; sp < n_spin; ++sp){
      G_tau[sp][tau] = get_G_tau(double(tau));
    }
  }

  return G_tau;
}
 

inline std::tuple<solver_core, solve_params_t, u_tau_t> test_setup(int n_site, int n_bath, int n_spin, double U, double mu, double t,
                                                                   constr_params_t const &cp, mat_t const &theta, vec_t const &epsilon) {
  // Set up the Solver
  solver_core S(cp);

  // create hybridization:
  for (auto const &tau : S.Delta_tau[0].mesh()) {
    double val;

    for (int block = 0; block < cp.gf_struct.size(); block++) {
      S.Delta_tau[block][tau] = 0.0;
      for (int i = 0; i < n_site; i++) {
        for (int j = 0; j < n_site; j++) {
          for (int n = 0; n < n_bath; n++) {
            S.Delta_tau[block][tau](i, j) += theta(i, n) * dagger(theta)(n, j) * one_fermion(tau, epsilon(n), cp.beta);
          }
        }
      }
    }
  }

  for (int block = 0; block < cp.gf_struct.size(); block++)
    for (int i = 0; i < n_site; i++)
      for (int j = 0; j < n_site; j++) std::printf("%d %d % 4.8f\n", i, j, S.Delta_tau[block][cp.n_tau - 1](i, j));
  //exit(0);

  // h_imp initialization: Hamiltonian of the impurity sites (n_site)

  many_body_operator h_imp, h_hyb, h_bath;
  for (int j = 0; j < n_site; j++) {
    h_imp -= mu * n("up", j);

    if (n_spin == 2) {
      h_imp -= mu * n("dn", j);
      h_imp += U * n("up", j) * n("dn", j);
    }
    for (int i = 0; i < n_site; i++) {
      if (i != j) {
        h_imp -= t * c_dag("up", i) * c("up", j);
        if (n_spin == 2) h_imp -= t * c_dag("dn", i) * c("dn", j);
      }
    }
  }

  // fundamental operator sets initialization
  // bath + imp = tot
  auto [fops_tot, qn_tot]   = make_fops(n_site, n_bath, n_site, n_spin);
  auto [fops_imp, qn_imp]   = make_fops(n_site, 0, n_site, n_spin);
  auto [fops_bath, qn_bath] = make_fops(0, n_bath, n_site, n_spin);

  // Solve Parameters
  solve_params_t sp;
  sp.h_imp           = h_imp;
  sp.n_cycles        = 50000;
  sp.length_cycle    = 10;
  sp.n_warmup_cycles = 20;
  sp.max_time        = -1;
  sp.verbosity       = 3;
  sp.post_process    = true;
  sp.measure_sign    = true;
  sp.quantum_numbers = qn_imp;
  sp.random_seed     = 12345789 + 928374 * mpi::communicator().rank();

  //auto h_hyb  = 0.0 * n("up", 0);
  //auto h_bath = 0.0 * n("up", n_site);

  // define h_hyb alone
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

  // define h_bath alone
  for (int k = 0; k < n_bath; k++) {
    h_bath += epsilon(k) * n("up", k + n_site);
    if (n_spin == 2) h_bath += epsilon(k) * n("dn", k + n_site);
  }

  // define the 3 different atom_diag (ED calculation with Triqs):
  auto ad_tot  = triqs::atom_diag::atom_diag<false>(h_imp + h_bath + h_hyb, fops_tot);
  auto ad_imp  = triqs::atom_diag::atom_diag<false>(h_imp, fops_imp, qn_imp);
  auto ad_bath = triqs::atom_diag::atom_diag<false>(h_bath, fops_bath);

  u_tau_t u_tau = make_ED_propagator(ad_tot, ad_imp, ad_bath, cp.beta, cp.n_tau);

  return {S, sp, u_tau};
}

inline void solve_cthyb(solver_core S, solve_params_t const &sp, u_tau_t const &u_tau, double tau_max) {

  // Solve the impurity model using cthyb
  auto result_cthyb = S.solve_cthyb(sp, tau_max);
  auto const &cp    = S.constr_params;

  for (int bl = 0; bl < result_cthyb.frame.size(); bl++)
    EXPECT_ARRAY_NEAR(((matrix_t)u_tau[bl][cp.n_tau - 1]), result_cthyb.frame[bl],
                      0.05 * u_tau[0][cp.n_tau - 1](0, 0)); // U[0](0,0) is essentially always the biggest value
}

inline void solve_selfconsistent(solver_core S, solve_params_t const &sp, u_tau_t const &u_tau, double tau_split, double tau_max) {

  std::printf("\n##################\nexact U(beta):\n");
  print(u_tau, tau_max);

  // Solve the impurity model using the self-consistency approach
  auto result_sc = S.solve_self_consistently(sp, u_tau, tau_split, tau_max);
  auto const &cp = S.constr_params;

  for (int bl = 0; bl < result_sc.frame.size(); bl++)
    EXPECT_ARRAY_NEAR(((matrix_t)u_tau[bl][cp.n_tau - 1]), result_sc.frame[bl], 0.05 * u_tau[0][cp.n_tau - 1](0, 0));
}

void solve_green(solver_core S, solve_params_t const &sp, u_tau_t const &u_tau) {
  S.solve_green(sp, u_tau);
}
