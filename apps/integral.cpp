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

std::pair<int, int> find_index(const std::vector<int> &block_shape, int iota) {
  int running_sum    = 0;
  int subspace_index = 0;
  int iota_p1        = iota + 1;

  for (int bl = 0; bl < block_shape.size(); ++bl) {
    int new_sum = running_sum + block_shape[bl];

    // If newSum exceeds a, we know that the number is in the subspace represented by block_shape[i]
    if (new_sum >= iota_p1) {
      subspace_index = iota_p1 - running_sum - 1; // -1 to make it zero-indexed
      return std::make_pair(bl, subspace_index);
    }

    running_sum = new_sum;
  }

  return std::make_pair(-1, -1); // Return a pair of -1s if not found
}

class BuildConfig {
  public:
  BuildConfig(frame_t const &frame_zeroth_order, double tau_split, double tau_max, std::vector<std::vector<fop_t>> const &all_d_ops,
              std::vector<std::vector<fop_t>> const &all_d_dag_ops, std::vector<int> const &block_shape, constr_params_t const &cp,
              hyb_tau_t const &Delta_tau, atom_diag const &ad_imp, interpolator_t<scalar_t> const &u_interpolator)
     : frame_zeroth_order(frame_zeroth_order),
       tau_split(tau_split),
       tau_max(tau_max),
       all_d_ops(all_d_ops),
       all_d_dag_ops(all_d_dag_ops),
       block_shape(block_shape),
       cp(cp),
       Delta_tau(Delta_tau),
       ad_imp(ad_imp),
       u_interpolator(u_interpolator) {}

  config_t operator()(auto const &tau_d_list, auto const &tau_d_dag_list, auto const &iota_d_list, auto const &iota_d_dag_list) {
    auto config = config_t(frame_zeroth_order, cp.gf_struct, {0.0, tau_split});
    for (auto i : range(tau_d_list.size())) {
      auto [bl, subspace_index]         = find_index(block_shape, iota_d_list[i]);
      auto d                            = all_d_ops[bl][subspace_index];
      d.tau                             = tau_d_list[i];
      auto [bl_dag, subspace_index_dag] = find_index(block_shape, iota_d_dag_list[i]);
      auto d_dag                        = all_d_dag_ops[bl_dag][subspace_index_dag];
      d_dag.tau                         = tau_d_dag_list[i];
      config.d_bl_list[bl].push_back(d);
      config.d_dag_bl_list[bl_dag].push_back(d_dag);
      config.d_list.push_back(d);
      config.d_dag_list.push_back(d_dag);
      config.split_times.push_back(d.tau);
      config.split_times.push_back(d_dag.tau);
    }
    auto diagram = diagram::time_diagram_t{config, {tau_split}};
    print_configuration(diagram);
    auto hyb_mat   = diagram::hyb_matrix_t(diagram, Delta_tau);
    auto hyb_wight = inclusion_exclusion(diagram, hyb_mat);
    std::cout << "hyb_wight: " << hyb_wight << std::endl;
    auto frame      = make_frame(impurity_product(ad_imp, diagram, tau_max, tau_split, &u_interpolator)
                                 * impurity_product(ad_imp, diagram, tau_split, 0, &u_interpolator));
    auto imp_weight = norm(frame);
    std::cout << "imp_weight: " << imp_weight << std::endl;
    return config;
  }

  private:
  frame_t const &frame_zeroth_order;
  double tau_split;
  double tau_max;
  std::vector<std::vector<fop_t>> const &all_d_ops;
  std::vector<std::vector<fop_t>> const &all_d_dag_ops;
  std::vector<int> const &block_shape;
  constr_params_t const &cp;
  hyb_tau_t const &Delta_tau;
  atom_diag const &ad_imp;
  interpolator_t<scalar_t> const &u_interpolator;
};

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
  auto frame      = make_frame(impurity_product(ad_imp, diagram, tau_max, tau_split, &u_interpolator)
                               * impurity_product(ad_imp, diagram, tau_split, 0, &u_interpolator));
  auto imp_weight = norm(frame);
  std::cout << "imp_weight: " << imp_weight << std::endl;

  std::vector<int> block_shape = {};
  for (auto [bl, bl_pair] : enumerate(cp.gf_struct)) {
    auto [bl_name, bl_size] = bl_pair;
    block_shape.push_back(bl_size);
  }

  int order            = 2;
  auto tau_d_list      = std::vector<double>(order, 0.0);
  auto tau_d_dag_list  = std::vector<double>(order, 0.0);
  auto iota_d_list     = std::vector<int>(order, 0);
  auto iota_d_dag_list = std::vector<int>(order, 0);
  tau_d_list           = {tau_max * 0.2, tau_max * 0.5};
  tau_d_dag_list       = {tau_max * 0.9, tau_max * 0.95};
  iota_d_list          = {0, 0};
  iota_d_dag_list      = {0, 0};

  auto build_config =
     BuildConfig(frame_zeroth_order, tau_split, tau_max, all_d_ops, all_d_dag_ops, block_shape, cp, Delta_tau, ad_imp, u_interpolator);
  auto config_new = build_config(tau_d_list, tau_d_dag_list, iota_d_list, iota_d_dag_list);
  return 0;
}
