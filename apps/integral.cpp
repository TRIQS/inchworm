#include <iostream>
#include <vector>
#include <algorithm>

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

using namespace inchworm;

template <typename T> void print_vector(const std::vector<T> &vec) {
  for (auto v : vec) std::cout << v << ' ';
  std::cout << std::endl;
}

std::vector<std::pair<std::vector<int>, std::vector<int>>> get_all_order(std::vector<int> const &order_list) {
  std::vector<std::pair<std::vector<int>, std::vector<int>>> res{};
  int order = order_list.size() / 2;
  std::vector<int> indicator(2 * order);
  std::fill(indicator.begin(), indicator.begin() + order, 0);
  std::fill(indicator.begin() + order, indicator.end(), 1);
  do {
    std::vector<int> order_d_ls     = {};
    std::vector<int> order_d_dag_ls = {};
    for (int i = 0; i < 2 * order; ++i) {
      if (indicator[i] == 0) {
        order_d_ls.push_back(order_list[i]);
      } else {
        order_d_dag_ls.push_back(order_list[i]);
      }
    }
    // std::cout << "tau_d_ls size: " << tau_d_ls.size() << std::endl;
    // std::cout << "tau_d_dag_ls size: " << tau_d_dag_ls.size() << std::endl;
    // Display the generated sub-vectors
    std::cerr << "First: " << std::endl;
    print_vector(order_d_ls);
    std::cerr << "Second: " << std::endl;
    print_vector(order_d_dag_ls);
    std::cerr << "---\n";
    res.push_back(std::make_pair(order_d_ls, order_d_dag_ls));
  } while (std::next_permutation(indicator.begin(), indicator.end()));
  return res;
}

template <typename T> std::vector<T> getElements(const std::vector<int> &indices, const std::vector<T> &values) {
  std::vector<T> result;
  for (int index : indices) {
    if (index >= 0 && index < values.size()) {
      result.push_back(values[index]);
    } else {
      // Handle out-of-range indices according to your requirements
      // For now, simply skipping them
    }
  }
  return result;
}

std::vector<double> changeVariable(const std::vector<double> &nus, double tau_max) {
  std::vector<double> taus(nus.size());
  taus[0] = nus[0] * tau_max;

  for (size_t it = 1; it < taus.size(); ++it) { taus[it] = taus[it - 1] + nus[it] * (tau_max - taus[it - 1]); }

  return taus;
}

double jacobian(const std::vector<double> &taus, double tau_max) {
  double prod = tau_max;
  for (size_t j = 1; j < taus.size(); ++j) { prod *= (tau_max - taus[j - 1]); }
  return std::abs(prod);
}

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
  BuildConfig(frame_t frame_zeroth_order, double tau_split, double tau_max, std::vector<std::vector<fop_t>> const &all_d_ops,
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
       u_interpolator(u_interpolator),
       config(config_t(frame_zeroth_order, cp.gf_struct, {0.0, tau_split})) {}

  void operator()(auto const &tau_d_list, auto const &tau_d_dag_list, auto const &iota_d_list, auto const &iota_d_dag_list) {
    clear_config();
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
    // auto diagram = diagram::time_diagram_t{config, {tau_split}};
    // auto hyb_mat = diagram::hyb_matrix_t(diagram, Delta_tau);
    // sign = diagram.sign();
    // hyb_weight   = inclusion_exclusion(diagram, hyb_mat);
    // u_products   = make_frame(impurity_product(ad_imp, diagram, tau_max, tau_split, &u_interpolator)
    //                           * impurity_product(ad_imp, diagram, tau_split, 0, &u_interpolator));
    // for (auto u_products_bl : u_products) { print_matrix(u_products_bl); }
    // auto imp_weight = norm(u_products);
    // print_configuration(diagram);
    // std::cout << "hyb_weight: " << hyb_weight << std::endl;
    // std::cout << "imp_weight: " << imp_weight << std::endl;
  }

  void evaluate_hyb_weight() {
    auto diagram = diagram::time_diagram_t{config, {tau_split}};
    auto hyb_mat = diagram::hyb_matrix_t(diagram, Delta_tau);
    sign         = diagram.sign();
    hyb_weight   = inclusion_exclusion(diagram, hyb_mat);
  }

  void evaluate_u_products() {
    auto diagram = diagram::time_diagram_t{config, {tau_split}};
    u_products   = make_frame(impurity_product(ad_imp, diagram, tau_max, tau_split, &u_interpolator)
                              * impurity_product(ad_imp, diagram, tau_split, 0, &u_interpolator));
  }

  public:
  hyb_scalar_t hyb_weight;
  frame_t u_products;
  double sign;

  private:
  config_t config;
  frame_t frame_zeroth_order;
  double tau_split;
  double tau_max;
  std::vector<std::vector<fop_t>> const &all_d_ops;
  std::vector<std::vector<fop_t>> const &all_d_dag_ops;
  std::vector<int> const &block_shape;
  constr_params_t const &cp;
  hyb_tau_t const &Delta_tau;
  atom_diag const &ad_imp;
  interpolator_t<scalar_t> const &u_interpolator;

  void clear_config() { config = config_t(frame_zeroth_order, cp.gf_struct, {0.0, tau_split}); }
};

int main() {

  constr_params_t cp;
  cp.beta        = 2.0;
  cp.gf_struct   = {{"up", 1}};
  cp.n_tau_green = 5;
  cp.n_tau_inch  = 21;
  cp.n_tau       = 10001;

  mat_t theta   = {{1.0}};
  vec_t epsilon = {1.0};

  int n_site = 1;
  int n_bath = epsilon.size();
  int n_spin = cp.gf_struct.size();
  double U   = 0.0;
  double mu  = 0.0;
  double t   = 0.0;

  auto [Delta_tau, ad_imp, u_tau, G_tau] = test_setup(n_site, n_bath, n_spin, U, mu, t, cp, theta, epsilon);

  double tau_max   = cp.beta;
  double tau_split = cp.beta * 0.855;
  std::cout << "Delta_tau shape:" << std::endl;
  print_block_shape(Delta_tau);
  std::cout << "G_tau shape:" << std::endl;
  print_block_shape(G_tau);
  std::cout << "u_tau shape:" << std::endl;
  print_block_shape(u_tau);
  std::cout << "ad_imp.n_subspaces(): " << ad_imp.n_subspaces() << std::endl;
  auto u_interpolator        = interpolator_t<scalar_t>(u_tau, u_tau[0].mesh().size());
  frame_t frame_zeroth_order = u_interpolator(tau_max - tau_split) * u_interpolator(tau_split);

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
  // // manually calculate the weight
  // auto config = config_t(frame_zeroth_order, cp.gf_struct, {0.0, tau_split});
  // auto d      = all_d_ops[1][0];
  // auto d_dag  = all_d_dag_ops[1][1];
  // d.tau       = tau_max * 0.2;
  // d_dag.tau   = tau_max * 0.95;
  // config.try_insert(d_dag, d);
  // d     = all_d_ops[0][0];
  // d_dag = all_d_dag_ops[0][1];
  // d.tau      = tau_max * 0.4;
  // d_dag.tau  = tau_max * 0.8;
  // config.try_insert(d_dag, d);
  // d     = all_d_ops[0][1];
  // d_dag = all_d_dag_ops[0][1];
  // d.tau      = tau_max * 0.98;
  // d_dag.tau  = tau_max * 0.1;
  // config.try_insert(d_dag, d);

  // auto diagram = diagram::time_diagram_t{config, {tau_split}};
  // print_configuration(diagram);
  // auto hyb_mat   = diagram::hyb_matrix_t(diagram, Delta_tau);
  // auto hyb_wight = inclusion_exclusion(diagram, hyb_mat);
  // std::cout << "hyb_wight: " << hyb_wight << std::endl;
  // auto frame      = make_frame(impurity_product(ad_imp, diagram, tau_max, tau_split, &u_interpolator)
  //                              * impurity_product(ad_imp, diagram, tau_split, 0, &u_interpolator));
  // auto imp_weight = norm(frame);
  // std::cout << "imp_weight: " << imp_weight << std::endl;

  //generate configuration from function
  std::vector<int> block_shape = {};
  for (auto [bl, bl_pair] : enumerate(cp.gf_struct)) {
    auto [bl_name, bl_size] = bl_pair;
    block_shape.push_back(bl_size);
  }
  std::vector<int> iota_d_list       = {0,0,0};
  std::vector<int> iota_d_dag_list   = {0,0,0};

  auto build_config =
     BuildConfig(frame_zeroth_order, tau_split, tau_max, all_d_ops, all_d_dag_ops, block_shape, cp, Delta_tau, ad_imp, u_interpolator);

  auto get_hyb_weight_sign = [&build_config, &iota_d_list, &iota_d_dag_list](auto const &tau_d_list, auto const &tau_d_dag_list) {
    build_config(tau_d_list, tau_d_dag_list, iota_d_list, iota_d_dag_list);
    build_config.evaluate_hyb_weight();
    return build_config.hyb_weight * build_config.sign;
  };

  auto get_u_tau_max_00 = [&build_config, &iota_d_list, &iota_d_dag_list](auto const &tau_d_list, auto const &tau_d_dag_list) {
    auto my_config = build_config;
    my_config(tau_d_list, tau_d_dag_list, iota_d_list, iota_d_dag_list);
    my_config.evaluate_hyb_weight();
    my_config.evaluate_u_products();
    if (my_config.u_products[0].size() == 0) return 0.0;
    return trace(my_config.u_products)* my_config.hyb_weight * my_config.sign;
  };
  // test decomposition of u_products_00
  auto [nui, wi]          = xfac::grid::QuadratureGK15(0, 1);
  int bondDim             = 40;
  int sweepBound          = 10;
  int n                   = 6;
  double integral         = 0.0;
  std::vector<int> pivot1 = {0,0,0, 14, 14,14};
  
  std::vector<double> nu1;
  for (int i = 0; i < pivot1.size(); i++) { nu1.push_back(nui[pivot1[i]]); }
  std::cout << "nu1: ";
  print_vector(nu1);
  std::vector<int> range(n);
  std::iota(range.begin(), range.end(), 0);
  auto order_list_pair = get_all_order(range);
  long count           = 0;
  auto f               = [&get_u_tau_max_00, &tau_max, &count, &order_list_pair](std::vector<double> nus) {
    auto taus  = changeVariable(nus, tau_max);
    double sum = 0.0;
    for (auto [order_c_list, order_c_dag_list] : order_list_pair) {
      sum += get_u_tau_max_00(getElements(order_c_list, taus), getElements(order_c_dag_list, taus));
    }
    count++;
    double j = jacobian(taus, tau_max);
    // std::cerr << "end taus:" << std::endl;
    // print_vector(taus);
    // std::cerr <<std::endl;
    return sum * j;
  };
  std::cout << std::setprecision(17) << "f(nu1): " << f(nu1) << std::endl;
  std::cout << "tci1" << std::endl;
  auto ci1 = xfac::CTensorCI<double, double>(f, std::vector(n, nui), {.pivot1 = pivot1});
  for (int i = 0; i < sweepBound; i++) {
    ci1.iterate();
    integral = ci1.sumWeighted(std::vector(n, wi));
    std::cout << i << " " << count << " " << ci1.pivotError[ci1.pivotError.size() - 1] << " " << integral << std::endl;
  }
  // std::cout << "tci2" << std::endl;
  // auto ci2 = xfac::CTensorCI2<double, double>(f, std::vector(n, nui), {.bond_dim = bondDim, .pivot1 = pivot1});
  // for (int i = 0; i < sweepBound; i++) {
  //   ci2.iterate();
  //   if (i == sweepBound - 1) { ci2.makeCanonical(); }
  //   integral = ci2.tt.sum(std::vector(n, wi));
  //   std::cout << i << " " << count << " " << ci2.pivotError[ci2.pivotError.size() - 1] << " " << integral << std::endl;
  // }
  std::cout << "frame_order_zeroth 00: "<< frame_zeroth_order[0](0,0) << std::endl;
  std::cout << "zero_order + first order "<< frame_zeroth_order[0](0,0) + integral<< std::endl;
  std::cout << "u_tau_max_00: "<< u_interpolator(tau_max)[0](0,0) << std::endl;

  return 0;
}
