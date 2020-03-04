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
    for (int j = 0; j < second_dim(m); j++) { std::printf("% 5.10f ", m(i, j)); }
    std::printf("\n");
  }
  std::printf("\n\n");
}

// Prepare funcdamental operator set
fundamental_operator_set make_fops(int N) {
  fundamental_operator_set fops;
  for (int o : range(N)) { fops.insert("up", o); }
  return fops;
}

TEST(inchworm, HubbardAtom) { // NOLINT

  // System Parameters
  double mu = 0.;
  //double h  = 0.1;

  // Construct Parameters
  constr_params_t cp;
  cp.beta      = 2.0;
  int n_site   = 2;
  cp.gf_struct = {{"up", {0}}, {"up", {1}}};
  cp.n_tau     = 500;
  cp.n_iw      = 250;

  // Set up the Solver
  solver_core S(cp);
  //int up = 0, dn = 1;
  int n_bath       = 2;
  double theta[]   = {3.0, 3.0};
  double epsilon[] = {0.0, 0.0};
  for (auto const &tau : S.Delta_tau[0].mesh()) {
    for (int i = 0; i < n_site; i++) {
      S.Delta_tau[i][tau] = 0.0;
      //for (int n = 0; n < n_bath; n++) {
      double val;
      if (epsilon[i] > 0.0)
        val = -theta[i] * theta[i] * (std::exp(-((double)tau) * (epsilon[i])) / (1. + std::exp(-cp.beta * epsilon[i])));
      else
        val = -theta[i] * theta[i] * (std::exp(-((double)tau - cp.beta) * (epsilon[i])) / (1. + std::exp(cp.beta * epsilon[i])));
      S.Delta_tau[i][tau] += val;
      //std::printf("hyb = %f, %f \n", val, (double)tau);
      //}
    }
  }

  std::cout << S.Delta_tau[0];

  std::vector<many_body_op_t> qn;
  qn.resize(1);
  auto h_int = 0 * n("up", 0);
  for (int j = 0; j < n_site; j++) {
    qn[0] += n("up", j);
    h_int -= mu * n("up", j);
  }

  // Solve Parameters
  solve_params_t sp;
  sp.h_int           = h_int;
  sp.n_cycles        = 500000;
  sp.length_cycle    = 10;
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
  auto fops      = make_fops(n_bath + n_site);
  auto fops_bath = make_fops(n_bath);

  auto h      = 0 * n("up", 0);
  auto h_bath = 0 * n("up", 0);
  for (int j = 0; j < n_site; j++) {
    h -= mu * (n("up", j) + n("dn", j));

    //for (int i = 0; i < n_bath; i++) {
      h += theta[j] * (c_dag("up", j) * c("up", j + n_site) + c_dag("up", j + n_site) * c("up", j));
      h += epsilon[j] * n("up", j + n_site);

      h_bath += epsilon[j] * n("up", j);
    //}
  }

  std::printf("\n\n");
  auto dtau = cp.beta;
  auto ad   = triqs::atom_diag::atom_diag<false>(h, fops);
  auto E0   = ad.get_gs_energy();
  auto ps   = partial_sum(ad, n_site, [dtau, E0](double E) { return std::exp(-dtau * (E - E0)); });
  //print_matrix(ps);

  auto ad_bath = triqs::atom_diag::atom_diag<false>(h_bath, fops_bath);
  auto ps_bath = partial_sum(ad_bath, 0, [dtau, E0](double E) { return std::exp(-dtau * (E - E0)); });
  //print_eigensystems(ad_bath);
  print_matrix(ps_bath);
  print_matrix(ps / ps_bath(0, 0));

  double x = theta[0] * cp.beta;
  /*  double tmp         = std::cosh(theta[0] * cp.beta / sqrt(2.));
  double total_serie = std::pow(tmp, 2);

  double order1 = (1. / 2.) * std::pow(x, 2);
  double order2 = (1. / 12.) * std::pow(x, 4);
  double order3 = (1. / 180.) * std::pow(x, 6);
  double order4 = (1. / 5040.) * std::pow(x, 8);
  */

  double tmp         = std::cosh(theta[0] * cp.beta / 2.);
  double total_serie = std::pow(tmp, 4);

  double order1 = (1. / 2.) * std::pow(x, 2);
  double order2 = (5. / 48.) * std::pow(x, 4);
  double order3 = (17. / 1440.) * std::pow(x, 6);
  double order4 = (13. / 16128.) * std::pow(x, 8);

  //beta**10*theta**10/7257600 + beta**8*theta**8/80640 + beta**6*theta**6/1440 + beta**4*theta**4/48 + beta**2*theta**2/4 + 1

  //std::printf("order 0: % 4.8f \n", order0 / zb);
  std::printf("order 1: % 4.8f \n", order1);
  std::printf("order 2: % 4.8f \n", order2);
  std::printf("order 3: % 4.8f \n", order3);
  std::printf("order 4: % 4.8f \n", order4);
  std::printf("\n\ntotal : % 4.8f \n", total_serie);
  //std::printf("\ntotal: % 4.8f \n", U00 / zb);
  //std::printf("\ntotal: % 4.8f \n", U11 / zb);
}

MAKE_MAIN
