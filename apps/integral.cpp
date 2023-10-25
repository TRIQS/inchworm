#include <iostream>
#include <vector>

#include <xfac/grid.h>
#include <xfac/tensor/tensor_ci.h>
#include <xfac/tensor/tensor_ci_2.h>
#include <xfac/tensor/tensor_train.h>
#include <inchworm/diagram/diagram.hpp>
#include <inchworm/diagram/inclusion_exclusion.hpp>
#include <inchworm/diagram/print.hpp>
#include <inchworm/atom_diag.hpp>
#include <inchworm/u_frame.hpp>
#include <inchworm/impurity_product.hpp>
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
  double tau_split = cp.beta * 0.8;
  std::cout << "Delta_tau shape:" << std::endl;
  print_block_shape(Delta_tau);
  std::cout << "G_tau shape:" << std::endl;
  print_block_shape(G_tau);
  std::cout << "u_tau shape:" << std::endl;
  print_block_shape(u_tau);
  std::cout << "ad_imp.n_subspaces(): " << ad_imp.n_subspaces() << std::endl;

  auto u_interpolator = interpolator_t<scalar_t>(u_tau, u_tau[0].mesh().size());

  frame_t frame_zeroth_order = u_interpolator(tau_max - tau_split) * u_interpolator(tau_split);

  auto config = config_t(frame_zeroth_order, cp.gf_struct, {0.0, tau_split});

  long n_bl                                     = cp.gf_struct.size();
  std::vector<std::vector<fop_t>> all_d_ops     = {};
  std::vector<std::vector<fop_t>> all_d_dag_ops = {};
  all_d_ops.resize(n_bl, {});
  all_d_dag_ops.resize(n_bl, {});
  auto fops = fundamental_operator_set{cp.gf_struct};

  for (auto [bl, bl_pair] : enumerate(cp.gf_struct)) {
    auto [bl_name, bl_size] = bl_pair;
    all_d_ops[bl].clear();
    all_d_dag_ops[bl].clear();
    for (auto idx : range(bl_size)) {
      all_d_ops[bl].push_back({0.0, false, fops[{bl_name, idx}], bl, idx});
      all_d_dag_ops[bl].push_back({0.0, true, fops[{bl_name, idx}], bl, idx});
    }
  }

  auto d     = all_d_ops[0][0];
  auto d_dag = all_d_dag_ops[0][0];
  d.tau      = tau_max * 0.5;
  d_dag.tau  = tau_max * 0.9;
  config.try_insert(d_dag, d);
  d.tau     = tau_max * 0.2;
  d_dag.tau = tau_max * 0.95;
  config.try_insert(d_dag, d);

  auto diagram = diagram::time_diagram_t{config, {tau_split}};
  print_configuration(diagram);
  auto hyb_mat   = diagram::hyb_matrix_t(diagram, Delta_tau);
  auto hyb_wight = inclusion_exclusion(diagram, hyb_mat);
  std::cout << "hyb_wight: " << hyb_wight << std::endl;
  auto frame = make_frame(impurity_product(ad_imp, diagram, tau_max, tau_split, &u_interpolator) * impurity_product(ad_imp, diagram, tau_split, 0, &u_interpolator));
  auto imp_weight = norm(frame);
  std::cout << "imp_weight: " << imp_weight << std::endl;

  return 0;
}
