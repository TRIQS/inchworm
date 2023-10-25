#include <iostream>
#include <vector>

#include <xfac/grid.h>
#include <xfac/tensor/tensor_ci.h>
#include <xfac/tensor/tensor_ci_2.h>
#include <xfac/tensor/tensor_train.h>
#include <inchworm/diagram/diagram.hpp>
#include <inchworm/atom_diag.hpp>
#include <inchworm/u_frame.hpp>
#include <inchworm/util.hpp>
#include <inchworm/interpolator.hpp>
#include "./hubbard.hpp"

using namespace xfac;
using namespace inchworm;

template <typename T> void print_block_shape(block_gf<imtime, T> const &x_tau) {
  std::cout << "number of taus: " << x_tau[0].mesh().size() << std::endl;
  std::cout << "number of block: " << x_tau.size() << std::endl;
  for (int bl = 0; bl < x_tau.size(); ++bl) { std::cout << "block: " << bl << " shape: " << x_tau[bl].target_shape() << std::endl; }
}

int main() {

  constr_params_t cp;
  cp.beta        = 2.0;
  cp.gf_struct   = {{"up", 2}, {"dn", 2}};
  cp.n_tau_green = 5;
  cp.n_tau_inch  = 21;
  cp.n_tau       = 10001;

  mat_t theta   = {{0.1, 0.3, 0.4}, {0.1, 0.2, 0.4}};
  vec_t epsilon = {1.0, -1.0, 1.2};

  int n_site = 2;
  int n_bath = epsilon.size();
  int n_spin = cp.gf_struct.size();
  double U   = 1.0;
  double mu  = 1.0;
  double t   = 1.0;

  auto [Delta_tau, ad_imp, u_tau, G_tau] = test_setup(n_site, n_bath, n_spin, U, mu, t, cp, theta, epsilon);

  double tau_max   = cp.beta;
  double tau_split = cp.beta * 0.9;
  std::cout << "Delta_tau shape:" << std::endl;
  print_block_shape(Delta_tau);
  std::cout << "G_tau shape:" << std::endl;
  print_block_shape(G_tau);
  std::cout << "u_tau shape:" << std::endl;
  print_block_shape(u_tau);
  std::cout << "ad_imp.n_subspaces(): " << ad_imp.n_subspaces() << std::endl;

  auto u_interpolator = interpolator_t<scalar_t>(u_tau, u_tau[0].mesh().size());

  frame_t frame_zeroth_order = u_interpolator(tau_max - tau_split) * u_interpolator(tau_split);


  return 0;
}
