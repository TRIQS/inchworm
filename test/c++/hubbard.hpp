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
#include <h5/h5.hpp>
#include <triqs/atom_diag/gf.hpp>

#include <triqs/test_tools/gfs.hpp>

using namespace inchworm;
using mat_t = triqs::arrays::matrix<double>;
using vec_t = triqs::arrays::array<double, 1>;

// Prepare fundamental operator set
inline std::pair<fundamental_operator_set, std::vector<many_body_op_t>> make_fops(int n_site, int n_bath, int bath_offset, int n_spin) {
  fundamental_operator_set fops;
  std::vector<many_body_op_t> qn(1);
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

// Exact non-interacting Green function of a single atomic level
double one_fermion(double tau, double eps, double beta) {
  if (eps >= 0) {
    return -std::exp(-tau * eps) / (1. + std::exp(-beta * eps));
  } else {
    return -std::exp((beta - tau) * eps) / (1. + std::exp(beta * eps));
  }
}

inline std::tuple<solver_core, solve_params_t, u_tau_t, g_tau_t> test_setup(int n_site, int n_bath, int n_spin, double U, double mu, double t,
                                                                            constr_params_t const &cp, mat_t const &theta, vec_t const &eps) {

  // === Define fundamental operator sets

  // Full system
  auto [fops_tot, qn_tot]   = make_fops(n_site, n_bath, n_site, n_spin);

  // Impurity
  auto [fops_imp, qn_imp]   = make_fops(n_site, 0, n_site, n_spin);

  // Bath
  auto [fops_bath, qn_bath] = make_fops(0, n_bath, n_site, n_spin);


  // === Initialize Hamiltonians

  // h_imp: Hamiltonian of the impurity sites (n_site)
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

  // h_bath: Hamiltonian of the bath (n_site)
  for (int k = 0; k < n_bath; k++) {
    h_bath += eps(k) * n("up", k + n_site);
    if (n_spin == 2) h_bath += eps(k) * n("dn", k + n_site);
  }

  // h_hyb: Hamiltonian coupling the impurity and the bath
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


  // === Set up the Solver

  solver_core S(cp);

  many_body_operator hyb_effective;
  
  // create hybridization:
  for (auto const &tau : S.Delta_tau[0].mesh()) {
    for (int block = 0; block < cp.gf_struct.size(); block++) {
      S.Delta_tau[block][tau] = 0.0;
      for (auto [i, j, n] : product_range(n_site, n_site, n_bath)) {
        S.Delta_tau[block][tau](i, j) += theta(i, n) * dagger(theta)(n, j) * one_fermion(tau, eps(n), cp.beta);
        
        if(tau==0){
          hyb_effective += 1. * c_dag("up", i) * c("up", j) ;
        }
      }
    }
  }


  // === Define the 3 different atom_diag objects (ED calculation with Triqs)

  auto ad_tot  = triqs::atom_diag::atom_diag<false>(h_imp + h_bath + h_hyb, fops_tot);
  //auto ad_imp  = triqs::atom_diag::atom_diag<false>(h_imp, fops_imp, qn_imp);
  auto ad_imp  = triqs::atom_diag::atom_diag<false>(h_imp, hyb_effective, fops_imp);
  auto ad_bath = triqs::atom_diag::atom_diag<false>(h_bath, fops_bath);

  // Calculate exact propagator
  u_tau_t u_tau = make_ED_propagator(ad_tot, ad_imp, ad_bath, cp.beta, cp.n_tau_inch);

  // Calculate exact Green function
  g_tau_t g_tau = real(atomic_g_tau(ad_tot, cp.beta, cp.gf_struct, cp.n_tau_green));



  // Solve Parameters
  solve_params_t sp;
  sp.h_imp           = h_imp;
  sp.n_cycles        = 500000;
  sp.length_cycle    = 10;
  sp.n_warmup_cycles = 20;
  sp.hyb_effective   = hyb_effective;
  sp.verbosity       = 10;


  return {S, sp, u_tau, g_tau};
}
