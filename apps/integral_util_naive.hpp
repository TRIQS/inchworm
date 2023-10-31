#include <boost/math/quadrature/gauss_kronrod.hpp>

template <int n> inline auto QuadratureGK(double a = 0, double b = 1) {
  static const auto abscissa = boost::math::quadrature::gauss_kronrod<double, n>::abscissa();
  static const auto weights  = boost::math::quadrature::gauss_kronrod<double, n>::weights();
  int nq                     = 2 * abscissa.size() - 1;
  std::vector<double> xi(nq), weight(nq);
  double factor = 0.5 * (b - a);
  for (uint i = 0; i < abscissa.size(); i++) {
    xi[nq / 2 + i]     = factor * (abscissa[i] + 1) + a;
    weight[nq / 2 + i] = weight[nq / 2 - i] = weights[i] * factor;
    xi[nq / 2 - i]                          = factor * (-abscissa[i] + 1) + a;
  }
  return make_pair(xi, weight);
}

template <typename T> void print_vector(const std::vector<T> &vec) {
  for (auto v : vec) std::cout << v << ' ';
  std::cout << std::endl;
}

inline std::vector<std::pair<std::vector<int>, std::vector<int>>> get_all_order(std::vector<int> const &order_list) {
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
    // std::cerr << "First: " << std::endl;
    // print_vector(order_d_ls);
    // std::cerr << "Second: " << std::endl;
    // print_vector(order_d_dag_ls);
    // std::cerr << "---\n";
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

inline std::vector<double> changeVariable(const std::vector<double> &nus, double tau_max, double tau_min = 0.0) {
  std::vector<double> taus(nus.size());
  taus[0] = nus[0] * (tau_max - tau_min) + tau_min;

  for (size_t it = 1; it < taus.size(); ++it) { taus[it] = taus[it - 1] + nus[it] * (tau_max - taus[it - 1]); }

  return taus;
}

inline double jacobian(const std::vector<double> &taus, double tau_max, double tau_min = 0.0) {
  double prod = tau_max - tau_min;
  for (size_t j = 1; j < taus.size(); ++j) { prod *= (tau_max - taus[j - 1]); }
  return std::abs(prod);
}

template <typename T> void print_block_shape(block_gf<imtime, T> const &x_tau) {
  std::cout << "number of taus: " << x_tau[0].mesh().size() << std::endl;
  std::cout << "number of block: " << x_tau.size() << std::endl;
  for (int bl = 0; bl < x_tau.size(); ++bl) { std::cout << "block: " << bl << " shape: " << x_tau[bl].target_shape() << std::endl; }
}

inline std::pair<int, int> find_index(const std::vector<int> &block_shape, int iota) {
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
  BuildConfig(frame_t &frame_zeroth_order, double tau_split, double tau_max, std::vector<std::vector<fop_t>> const &all_d_ops,
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
  }

  void evaluate_hyb_weight() {
    auto diagram = diagram::time_diagram_t{config, {tau_split}}; //FIXME: calculate diagram twice, diagram does not has default constructor
    auto hyb_mat = diagram::hyb_matrix_t(diagram, Delta_tau);
    sign         = diagram.sign();
    hyb_weight   = inclusion_exclusion(diagram, hyb_mat);
  }

  void evaluate_u_products() {
    auto diagram = diagram::time_diagram_t{config, {tau_split}};
    u_products   = make_frame(impurity_product(ad_imp, diagram, tau_max, tau_split, &u_interpolator)
                              * impurity_product(ad_imp, diagram, tau_split, 0, &u_interpolator));
    if (u_products[0].size() != 0) {
      auto hyb_mat = diagram::hyb_matrix_t(diagram, Delta_tau);
      sign         = diagram.sign();
      hyb_weight   = inclusion_exclusion(diagram, hyb_mat);
    }
  }

  public:
  hyb_scalar_t hyb_weight;
  frame_t u_products;
  int sign;

  config_t config;
  frame_t &frame_zeroth_order;
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

double evaluate_u_tau_max_00(frame_t &frame_zeroth_order, double tau_split, double tau_max, std::vector<std::vector<fop_t>> const &all_d_ops,
                             std::vector<std::vector<fop_t>> const &all_d_dag_ops, std::vector<int> const &block_shape, constr_params_t const &cp,
                             hyb_tau_t const &Delta_tau, atom_diag const &ad_imp, interpolator_t<scalar_t> const &u_interpolator,
                             auto const &tau_d_list, auto const &tau_d_dag_list, auto const &iota_d_list, auto const &iota_d_dag_list) {
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
  auto diagram      = diagram::time_diagram_t{config, {tau_split}};
  auto u_products   = make_frame(impurity_product(ad_imp, diagram, tau_max, tau_split, &u_interpolator)
                                 * impurity_product(ad_imp, diagram, tau_split, 0, &u_interpolator));
  int sign          = 0;
  double hyb_weight = 0.0;
  if (u_products[0].size() != 0) {
    auto hyb_mat = diagram::hyb_matrix_t(diagram, Delta_tau);
    sign         = diagram.sign();
    hyb_weight   = inclusion_exclusion(diagram, hyb_mat);
    return u_products[0](0, 0) * hyb_weight * sign;
  } else {
    return 0.0;
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

//  auto build_config =
//      BuildConfig(frame_zeroth_order, tau_split, tau_max, all_d_ops, all_d_dag_ops, block_shape, cp, Delta_tau, ad_imp, u_interpolator);

//   auto get_hyb_weight_sign = [&build_config, &iota_d_list, &iota_d_dag_list](auto const &tau_d_list, auto const &tau_d_dag_list) {
//     build_config(tau_d_list, tau_d_dag_list, iota_d_list, iota_d_dag_list);
//     build_config.evaluate_hyb_weight();
//     return build_config.hyb_weight * build_config.sign;
//   };

//   auto get_u_tau_max_00 = [&build_config, &iota_d_list, &iota_d_dag_list](auto const &tau_d_list, auto const &tau_d_dag_list) {
//     auto my_config = build_config;
//     my_config(tau_d_list, tau_d_dag_list, iota_d_list, iota_d_dag_list);
//     my_config.evaluate_hyb_weight();
//     my_config.evaluate_u_products();
// std::cout << std::endl;
// std::cout << "split_times: " << std::endl;
// print_vector(my_config.config.split_times);
// std::cout << "sign: " << my_config.sign << std::endl;
// std::cout << "hyb_weight: " << my_config.hyb_weight << std::endl;
// std::cout << "u_products: " << std::endl;
// for (auto u_products_bl : my_config.u_products) { print_matrix(u_products_bl); }
//     if (my_config.u_products[0].size() == 0) return 0.0;
//     return my_config.u_products[0](0, 0) * my_config.hyb_weight * my_config.sign;
//   };