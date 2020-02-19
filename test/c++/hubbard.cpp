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
#include <triqs/atom_diag/partial_trace.hpp>
#include <inchworm/solver_core.hpp>

#include <triqs/gfs.hpp>
#include <triqs/h5.hpp>
#include <triqs/test_tools/gfs.hpp>

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
    for (int j = 0; j < second_dim(m); j++) { std::printf("% 5.8f ", m(i, j)); }
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

TEST(inchworm, HubbardAtom) { // NOLINT

  // System Parameters
  double U  = 0.;
  double mu = U / 2;
  //double h  = 0.1;

  // Construct Parameters
  constr_params_t cp;
  cp.beta      = 2.0;
  cp.gf_struct = {{"up", {0}}, {"dn", {0}}};
  cp.n_tau     = 3;
  cp.n_iw      = 1;

  // Set up the Solver
  solver_core S(cp);
  //int up = 0, dn = 1;
  int n_bath       = 1;
  double theta[]   = {0.1, 0.05, 0.05};
  double epsilon[] = {0.0, 0.0, 0.0};

  for (auto const &tau : S.Delta_tau[0].mesh()) {
    for (int i = 0; i < 2; i++) {
      S.Delta_tau[i][tau] = 0.0;
      for (int n = 0; n < n_bath; n++) {
        // factor 2 is for spin:
        S.Delta_tau[i][tau] -= 2 * theta[n] * theta[n] * (std::exp(-tau * epsilon[n]) / (1. + std::exp(-cp.beta * epsilon[n])));
      }
    }
  }

  std::cout << S.Delta_tau[0];

  std::vector<many_body_op_t> qn;
  qn.resize(1);
  qn[0] += n("up", 0) + n("dn", 0);

  // Solve Parameters
  solve_params_t sp;
  sp.h_int           = U * n("up", 0) * n("dn", 0) - mu * (n("up", 0) + n("dn", 0));
  sp.n_cycles        = 10000;
  sp.length_cycle    = 4;
  sp.n_warmup_cycles = 20;
  sp.max_time        = -1;
  sp.verbosity       = 3;
  sp.post_process    = true;
  sp.measure_sign    = true;
  sp.quantum_numbers = qn;

  // Solve the impurity model
  S.solve_single_step(sp);

  // Store the Result
  {
    auto arch = triqs::h5::file("hubbard.out.h5", 'w');
    h5_write(arch, "S", S);
  }

  // Compare against the reference data
  // h5diff("hubbard.out.h5", "hubbard.ref.h5")
  auto fops      = make_fops(n_bath + 1);
  auto fops_bath = make_fops(n_bath);

  auto h = U * (n("up", 0) * (n("dn", 0))); // 0 is the only interacting orbital
  h -= mu * (n("up", 0) + n("dn", 0));

  auto h_bath = 0 * (n("up", 0) + (n("dn", 0)));

  for (int i = 0; i < n_bath; i++) {
    h += theta[i] * (c_dag("up", 0) * c("up", i + 1) + c_dag("up", i + 1) * c("up", 0));
    h += theta[i] * (c_dag("dn", 0) * c("dn", i + 1) + c_dag("dn", i + 1) * c("dn", 0));
    h += epsilon[i] * (n("up", i + 1) + n("dn", i + 1));

    h_bath += epsilon[i] * (n("up", i) + n("dn", i));
    //h -= mu * (n("up", i + 1) + n("dn", i + 1));
  }

  //std::vector<many_body_op_t> qn;
  //qn.resize(1);
  //qn[0] += n("up", 0) + n("dn", 0);

  std::printf("\n\n");
  auto dtau = cp.beta;
  auto ad   = triqs::atom_diag::atom_diag<false>(h, fops);
  auto ps   = partial_sum(ad, 2, [dtau](double x) { return std::exp(-dtau * x); });
  print_matrix(ps);

  auto ad_bath = triqs::atom_diag::atom_diag<false>(h_bath, fops_bath);
  auto ps_bath = partial_sum(ad_bath, 0, [dtau](double x) { return std::exp(-dtau * x); });
  //print_eigensystems(ad_bath);
  //print_matrix(ps_bath);
  print_matrix(ps / ps_bath(0, 0));

  double analy1 = -(theta[0] * theta[0] * cp.beta * cp.beta);
  std::printf("\n\n  % 4.8f \n", analy1);

  double analy2 = (theta[0] * theta[0] * cp.beta * cp.beta * theta[0] * theta[0] * cp.beta * cp.beta);
  std::printf("\n\n  % 4.8f \n", analy2);
}

MAKE_MAIN
