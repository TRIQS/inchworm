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

// Prepare funcdamental operator set
fundamental_operator_set make_fops(int idx1, int idx2) {
  fundamental_operator_set fops;
  for (int i = idx1; i < idx2; i++) {
    fops.insert("up", i);
    //fops.insert("dn", i);
  }
  return fops;
}

TEST(inchworm, HubbardAtom) { // NOLINT

  // System Parameters
  double U  = 0.;
  double mu = U / 2.;
  //double h  = 0.1;

  // Construct Parameters
  constr_params_t cp;
  cp.beta      = 2.0;
  int n_site   = 1;
  cp.gf_struct = {{"up", {0}}};
  cp.n_tau     = 500;
  cp.n_iw      = 250;

  // Set up the Solver
  solver_core S(cp);
  //int up = 0, dn = 1;
  int n_bath       = 1;
  double theta[]   = {2.0, 0.5,0.5, 0.5,0.5, 0.5,0.5, 0.5};
  double epsilon[] = {4.0,-.00,.00,-.00,.00,-.00,.00,-.00};

  for (auto const &tau : S.Delta_tau[0].mesh()) {
    double val;
    for (int i = 0; i < 1 * n_site; i++) {
      S.Delta_tau[i][tau] = 0.0;
      for (int n = 0; n < n_bath; n++) {
        //if (epsilon[n] >= 0.0)
        //val = -theta[n] * theta[n] * (std::exp(-((double)tau) * (epsilon[n])) / (1. + std::exp(-cp.beta * epsilon[n])));
        //else
        val = -theta[n] * theta[n] * (std::exp(-((double)tau - cp.beta) * (epsilon[n])) / (1. + std::exp(cp.beta * epsilon[n])));
        S.Delta_tau[i][tau] += val;
      }
    }
  }

  std::vector<many_body_op_t> qn;
  qn.resize(1);
  auto h_atom = 0 * (n("up", 0));
  for (int j = 0; j < n_site; j++) {
    qn[0] += n("up", j);
    h_atom += -mu * (n("up", j));
  }

  // Solve Parameters
  solve_params_t sp;
  sp.h_int           = h_atom;
  sp.n_cycles        = 100000;
  sp.length_cycle    = 10;
  sp.n_warmup_cycles = 20;
  sp.max_time        = -1;
  sp.verbosity       = 3;
  sp.post_process    = true;
  sp.measure_sign    = true;
  sp.quantum_numbers = qn;
  sp.random_seed     = 22345789 + 928374 * mpi::communicator().rank();

  // Solve the impurity model
  S.solve_single_step(sp);
  //exit(0);

  // Compare against the reference data
  // h5diff("hubbard.out.h5", "hubbard.ref.h5")
  auto fops_tot  = make_fops(0, n_site + n_bath);
  auto fops_atom = make_fops(0, n_site);
  auto fops_bath = make_fops(n_site, n_site + n_bath);

  //std::printf("salut\n");
  auto h_hyb  = 0.0 * n("up", 0);
  auto h_bath = 0.0 * n("up", 0);
  for (int j = 0; j < n_site; j++) {
    for (int i = 0; i < n_bath; i++) {
      h_hyb += theta[i] * (c_dag("up", j) * c("up", i + n_site) + c_dag("up", i + n_site) * c("up", j));
      //h_hyb += theta[i] * (c_dag("dn", j) * c("dn", i + n_site) + c_dag("dn", i + n_site) * c("dn", j));

      h_bath += epsilon[i] * (n("up", i + n_site));
    }
  }

  //std::printf("\n\n");
  auto dtau    = cp.beta;
  auto ad_tot  = triqs::atom_diag::atom_diag<false>(h_atom + h_bath + h_hyb, fops_tot);
  auto ad_atom = triqs::atom_diag::atom_diag<false>(h_atom, fops_atom, qn);
  auto ad_bath = triqs::atom_diag::atom_diag<false>(h_bath, fops_bath);

  std::printf("\nad_tot:\n");
  //print_eigensystems(ad_tot);
  std::printf("\nad_atom:\n");
  //print_eigensystems(ad_atom);
  std::printf("\nad_bath:\n");
  //print_eigensystems(ad_bath);

  u_tau_t u_tau = make_ED_propagator(ad_tot, ad_atom, ad_bath, cp.beta, cp.n_tau);

  /*
  auto E0      = ad_tot.get_gs_energy();
  auto u_frame = partial_trace(ad_tot, ad_atom, [dtau, E0](double E) { return std::exp(-dtau * (E - E0)); });
  std::printf("\n");
  auto Z_bath = trace(ad_bath, [dtau, E0](double E) { return std::exp(-dtau * (E - E0)); });
  print(u_frame, 1. / Z_bath);

  u_tau_t u_tau = make_propagator(ad_atom, 3);
  assign_u_frame_to_propagator(u_tau, u_frame, 1, 1. / Z_bath);

  std::printf("\n");
  for (int i_tau = 0; i_tau < cp.n_tau; i_tau++) { print(u_tau, i_tau); }
  std::printf("\n");
*/

  double tau_split = 0.95 * cp.beta;
  double tau_max   = 1.00 * cp.beta;

  S.solve_self_consistently(sp, u_tau, tau_split, tau_max);
  //std::printf("\n\n");
  //print(u_tau, tau_split);
  std::printf("\nexact U(beta):\n");
  print(u_tau, tau_max);

  fprint(u_tau, cp.n_tau);

  double order_0[2];
  double order_1[2];
  double order_2[2];

  int N_tau = 4000;
  for (int bl = 0; bl < 2; bl++) {
    double integral = 0.0;
    double dtau1    = tau_split / (N_tau - 1.);
    double dtau2    = (tau_max - tau_split) / (N_tau - 1.);
    for (int i_tau = 0; i_tau < N_tau; i_tau++) {
      for (int j_tau = 0; j_tau < N_tau; j_tau++) {
        double tau1 = (tau_split)*i_tau / (N_tau - 1.);
        double tau2 = (tau_max - tau_split) * j_tau / (N_tau - 1.) + tau_split;

        //if ((tau1 <= tau_split) and (tau_split <= tau2) and (tau2 < tau_max)) {
        double factor = 1.0;
        if (0 == i_tau) factor *= 0.5;
        if (N_tau - 1 == i_tau) factor *= 0.5;
        if (0 == j_tau) factor *= 0.5;
        if (N_tau - 1 == j_tau) factor *= 0.5;

        //std::printf("tau1 =% f, tau2 =% f  \n", tau1, tau2);
        //std::printf("% f  % f  % f  % f \n", tau1, tau_split - tau1, tau2 - tau_split, tau_max - tau2);
        //fflush(stdout);

        //std::printf("salut\n"); fflush(stdout);
        double u_tau1 = u_tau[bl](tau1 - 0.0)(0, 0);
        double u_tau2 = u_tau[1 - bl](tau_split - tau1)(0, 0);
        double u_tau3 = u_tau[1 - bl](tau2 - tau_split)(0, 0);
        double u_tau4 = u_tau[bl](tau_max - tau2)(0, 0);

        double w_hyb;
        if (bl == 0)
          w_hyb = S.Delta_tau[0](tau2 - tau1)(0, 0);
        else
          w_hyb = -S.Delta_tau[0](cp.beta + tau1 - tau2)(0, 0);

        if (bl == 0) factor = -factor;
        integral += factor * dtau1 * dtau2 * u_tau1 * u_tau2 * u_tau3 * u_tau4 * w_hyb;
        //}
      }
    }
    order_0[bl] = (double)(u_tau[bl](tau_max - tau_split)(0, 0) * u_tau[bl](tau_split)(0, 0));
    order_1[bl] = integral;

    //std::printf("\norder0 = %f  \n", order0);
    //std::printf("order1 = %f  \n", integral);
    //std::printf("\nsum = %f  \n", integral + order0);
  }

  ///*

  N_tau = 80;
  for (int bl = 0; bl < 2; bl++) {
    double integral = 0.0;
    double dtau1    = tau_split / (N_tau - 1.);
    double dtau2    = (tau_max - tau_split) / (N_tau - 1.);
    for (int i_tau = 0; i_tau < N_tau; i_tau++) {
      for (int j_tau = 0; j_tau < N_tau; j_tau++) {
        for (int k_tau = 0; k_tau < N_tau; k_tau++) {
          for (int l_tau = 0; l_tau < N_tau; l_tau++) {
            double tau1  = (tau_split)*i_tau / (N_tau - 1.);
            double tau1p = (tau_split)*j_tau / (N_tau - 1.);
            double tau2  = (tau_max - tau_split) * k_tau / (N_tau - 1.) + tau_split;
            double tau2p = (tau_max - tau_split) * l_tau / (N_tau - 1.) + tau_split;

            if ((tau1 < tau1p) and (tau2 < tau2p)) {
              double factor = 1.0;
              if (0 == i_tau) factor *= 0.5;
              if (0 == j_tau) factor *= 0.5;
              if (0 == k_tau) factor *= 0.5;
              if (0 == l_tau) factor *= 0.5;
              if (N_tau - 1 == i_tau) factor *= 0.5;
              if (N_tau - 1 == j_tau) factor *= 0.5;
              if (N_tau - 1 == k_tau) factor *= 0.5;
              if (N_tau - 1 == l_tau) factor *= 0.5;
              //if (tau_split == tau1) factor *= 0.5;
              //if (tau_split == tau2) factor *= 0.5;
              //if (0.0 == tau1) factor *= 0.5;
              //if (tau_max == tau2) factor *= 0.5;
              //std::printf("tau1 =% f, tau2 =% f  \n", tau1, tau2);
              //std::printf("% f  % f  % f  % f \n", tau1, tau_split - tau1, tau2 - tau_split, tau_max - tau2);
              //fflush(stdout);

              //std::printf("salut\n"); fflush(stdout);
              double u_tau1 = u_tau[bl](tau1 - 0.0)(0, 0);
              double u_tau2 = u_tau[1 - bl](tau1p - tau1)(0, 0);
              double u_tau3 = u_tau[bl](tau_split - tau1p)(0, 0);
              double u_tau4 = u_tau[bl](tau2 - tau_split)(0, 0);
              double u_tau5 = u_tau[1 - bl](tau2p - tau2)(0, 0);
              double u_tau6 = u_tau[bl](tau_max - tau2)(0, 0);

              double w_hyb;
              if (bl == 0)
                w_hyb = S.Delta_tau[0](tau2p - tau1)(0, 0) * S.Delta_tau[0](cp.beta + tau1p - tau2)(0, 0);
              else
                w_hyb = S.Delta_tau[0](tau2 - tau1p)(0, 0) * S.Delta_tau[0](cp.beta + tau1 - tau2p)(0, 0);

              integral += factor * dtau1 * dtau1 * dtau2 * dtau2 * u_tau1 * u_tau2 * u_tau3 * u_tau4 * u_tau5 * u_tau6 * w_hyb;
            }
          }
        }
      }
    }
    order_2[bl] = integral;
    //std::printf("order2 = %f  \n", integral);
  } //*/

  for (int bl = 0; bl < 2; bl++) {
    std::printf("\norder0 = %f", order_0[bl]);
    std::printf("\norder1 = %f ", order_1[bl]);
    std::printf("\norder2 = %f ", order_2[bl]);
    std::printf("\nsum = %f  \n\n", order_0[bl] + order_1[bl] + order_2[bl]);
  }

  //*
  double B   = tau_max * theta[0] / 2.;
  double t_s = tau_split * theta[0] / 2.;

  auto y1 = [](double Dt) {
    double s2 = std::sinh(2 * Dt);
    double c2 = std::cosh(2 * Dt);
    return 1. / 8 * Dt * (c2 + 2.) + 5. / 16 * s2;
  };

  auto y2 = [](double Dt) {
    double s3 = std::sinh(2 * Dt);
    double c3 = std::cosh(2 * Dt);
    return 1. / 64 * Dt * Dt * (c3 + 4.) + 15. / 128 * Dt * s3 + 3. / 32 * (c3 - 1.);
  };

  auto y3 = [](double Dt) {
    double s4 = std::sinh(2 * Dt);
    double c4 = std::cosh(2 * Dt);
    return 1. / 768 * Dt * Dt * Dt * (c4 + 8.) + 5. / 256 * Dt * Dt * s4 + 1. / 1024 * Dt * (57. * c4 - 64.) + 7. / 2048 * s4;
  };

  auto y4 = [](double Dt) {
    double s5 = std::sinh(2 * Dt);
    double c5 = std::cosh(2 * Dt);
    return 1. / 12288 * Dt * Dt * Dt * Dt * (c5 + 16.) + 25. / 12288 * Dt * Dt * Dt * s5 + 1. / 16384 * Dt * Dt * (205. * c5 - 320.)
       + 435. / 32768 * Dt * s5 - 5. / 512. * (c5 - 1.);
  };

  auto y5 = [](double Dt) {
    double s6 = std::sinh(2 * Dt);
    double c6 = std::cosh(2 * Dt);
    return 1. / 245760 * Dt * Dt * Dt * Dt * Dt * (c6 + 32.) + 5. / 32768 * Dt * Dt * Dt * Dt * s6 + 1. / 65536 * Dt * Dt * Dt * (107. * c6 - 256.)
       + 316. / 65536 * Dt * Dt * s6 - 1. / 262144 * Dt * (279. * c6 - 2304.) - 2025. / 524288 * s6;
  };

  double tmp    = std::cosh(B - t_s) * std::cosh(t_s - 0.0);
  double order0 = tmp * tmp;
  std::printf("order 0:       % 4.8f \n", order0);

  double order1 = 2 * y1(t_s) * y1(B - t_s);
  std::printf("order 1:       % 4.8f \n", order1);

  double order2 = 4 * (y2(t_s) * y2(B - t_s));
  std::printf("order 2:       % 4.8f \n", order2);

  double order3 = -8 * (y1(t_s) * y5(B - t_s) + y2(t_s) * y4(B - t_s) + y4(t_s) * y2(B - t_s) + y5(t_s) * y1(B - t_s));
  std::printf("order 3:       % 4.8f \n", order3);

  //std::printf("order 0+1+2+3: % 4.8f \n", order0 + order1*2. + order2/4.*6. + order3/36.*24);
  std::printf("order 0+1+2+3: % 4.8f \n", order0 + order1 + order2 + order3);
  std::printf("\n\ncosh^2(B): % 4.8f \n", std::cosh(B) * std::cosh(B));
  //*/
}

MAKE_MAIN
