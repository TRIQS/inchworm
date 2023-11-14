#pragma once
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

inline unsigned long long factorial(int n) {
  unsigned long long result = 1;
  for (int i = 1; i <= n; ++i) { result *= i; }
  return result;
}

template <int n> inline auto quadrature_GK(double a = 0, double b = 1) {
  static const auto abscissa = boost::math::quadrature::gauss_kronrod<double, n>::abscissa();
  static const auto weights  = boost::math::quadrature::gauss_kronrod<double, n>::weights();
  int nq                     = 2 * abscissa.size() - 1;
  std::vector<double> xi(nq), weight(nq);
  double factor = 0.5 * (b - a);
  for (unsigned int i = 0; i < abscissa.size(); i++) {
    xi[nq / 2 + i]     = factor * (abscissa[i] + 1) + a;
    weight[nq / 2 + i] = weight[nq / 2 - i] = weights[i] * factor;
    xi[nq / 2 - i]                          = factor * (-abscissa[i] + 1) + a;
  }
  return make_pair(xi, weight);
}

inline auto select_quadrature_GK(int n, double a = 0, double b = 1) {
  switch (n) {
    case 15: return quadrature_GK<15>(a, b);
    case 30: return quadrature_GK<30>(a, b);
    case 45: return quadrature_GK<45>(a, b);
    default: {
      std::cout << "Not supported\n";
      return quadrature_GK<15>(a, b);
    };
  }
}

template <typename T> void print_vector(const std::vector<T> &vec) {
  for (auto v : vec) std::cout << v << ' ';
  std::cout << std::endl;
}

inline std::vector<std::pair<std::vector<int>, std::vector<int>>> get_all_phi(std::vector<int> const &range) {
  std::vector<std::pair<std::vector<int>, std::vector<int>>> res{};
  int order = range.size() / 2;
  std::vector<int> indicator(2 * order);
  std::fill(indicator.begin(), indicator.begin() + order, 0);
  std::fill(indicator.begin() + order, indicator.end(), 1);
  do {
    std::vector<int> order_d_ls     = {};
    std::vector<int> order_d_dag_ls = {};
    for (int i = 0; i < 2 * order; ++i) {
      if (indicator[i] == 0) {
        order_d_ls.push_back(range[i]);
      } else {
        order_d_dag_ls.push_back(range[i]);
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

inline std::vector<std::pair<std::vector<int>, std::vector<int>>> get_all_phi_crossing(std::vector<int> const &range) {
  std::vector<std::pair<std::vector<int>, std::vector<int>>> res{};
  int order                       = range.size() / 2;
  std::vector<int> order_d_ls     = {};
  std::vector<int> order_d_dag_ls = {};
  for (int i = 0; i < 2 * order; i = i + 2) {
    order_d_ls.push_back(range[i]);
    order_d_dag_ls.push_back(range[i + 1]);
  }
  // std::cout << "tau_d_ls size: " << order_d_ls.size() << std::endl;
  // std::cout << "tau_d_dag_ls size: " << order_d_dag_ls.size() << std::endl;
  res.push_back(std::make_pair(order_d_ls, order_d_dag_ls));
  order_d_ls.clear();
  order_d_dag_ls.clear();
  for (int i = 0; i < 2 * order; i = i + 2) {
    order_d_ls.push_back(range[i + 1]);
    order_d_dag_ls.push_back(range[i]);
  }
  res.push_back(std::make_pair(order_d_ls, order_d_dag_ls));
  return res;
}

// geneerate all samples of B.size() from A and store them in all_combinations; each element of B can take any value from A
template <typename T>
void generate_combinations(const std::vector<T> &A, std::vector<T> &B, int idx, std::vector<std::vector<T>> &all_combinations) {
  if (idx == B.size()) {
    all_combinations.push_back(B);
    return;
  }

  for (int i = 0; i < A.size(); ++i) {
    B[idx] = A[i];
    generate_combinations(A, B, idx + 1, all_combinations);
  }
}

inline std::vector<int> generate_number_in_block(const std::vector<int> &block_shape, const std::vector<int> &iota_d_list) {
  std::vector<int> res(block_shape.size(), 0);
  for (auto iota : iota_d_list) {
    for (int bl = 0; bl < block_shape.size(); ++bl) {
      if (iota < block_shape[bl]) {
        res[bl] += 1;
        break;
      } else {
        iota -= block_shape[bl];
      }
    }
  }
  return res;
}

inline std::vector<std::pair<std::vector<int>, std::vector<int>>> get_all_iota(const std::vector<int> &block_shape, int order) {
  std::vector<std::pair<std::vector<int>, std::vector<int>>> res{};
  int n_phi = std::accumulate(block_shape.begin(), block_shape.end(), 0);
  std::vector<std::vector<int>> all_iota{};
  std::vector<int> iota(order, 0);
  std::vector<int> range(n_phi);
  std::iota(range.begin(), range.end(), 0);
  generate_combinations(range, iota, 0, all_iota);
  std::vector<std::vector<int>> all_number_in_block{};

  for (auto iota_d_list : all_iota) { all_number_in_block.push_back(generate_number_in_block(block_shape, iota_d_list)); }

  size_t i = 0;
  for (auto iota_d_list : all_iota) {
    size_t j = 0;
    for (auto iota_d_dag_list : all_iota) {
      if (all_number_in_block[i] == all_number_in_block[j]) { res.push_back(std::make_pair(iota_d_list, iota_d_dag_list)); }
      j++;
    }
    i++;
  }
  return res;
}

template <typename T> std::vector<T> get_elements(const std::vector<int> &indices, const std::vector<T> &values) {
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

// inline std::vector<double> change_variable(const std::vector<double> &nus, double tau_max, double tau_min = 0.0) {
//   std::vector<double> taus(nus.size());
//   taus[nus.size() - 1] = tau_min + (tau_max - tau_min) * std::pow(nus[nus.size() - 1], (1.0 / nus.size()));
//   for(int i = nus.size() - 2; i >= 0; --i) {
//     taus[i] = tau_min + (taus[i+1]- tau_min) * std::pow(nus[i], (1.0 /(i + 1)));
//   }
//   return taus;
// }

// inline double jacobian(const std::vector<double> &taus, double tau_max, double tau_min = 0.0) {
//    return std::pow(tau_max - tau_min, taus.size())/factorial(taus.size());
// }

inline std::vector<double> change_variable(const std::vector<double> &nus, double tau_max, double tau_min = 0.0) {
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

inline std::pair<int, int> findIndex(const std::vector<int> &block_shape, int iota) {
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

double evaluate_u_tau_max(frame_t &frame_zeroth_order, double tau_split, double tau_max, std::vector<std::vector<fop_t>> const &all_d_ops,
                          std::vector<std::vector<fop_t>> const &all_d_dag_ops, std::vector<int> const &block_shape, constr_params_t const &cp,
                          hyb_tau_t const &Delta_tau, atom_diag const &ad_imp, interpolator_t<scalar_t> const &u_interpolator, auto const &tau_d_list,
                          auto const &tau_d_dag_list, auto const &iota_d_list, auto const &iota_d_dag_list, int bl_indx, int subspace_indx) {
  auto config = config_t(frame_zeroth_order, cp.gf_struct, {0.0, tau_split});
  for (auto i : range(tau_d_list.size())) {
    auto [bl, subspace_d_index]       = findIndex(block_shape, iota_d_list[i]);
    auto d                            = all_d_ops[bl][subspace_d_index];
    d.tau                             = tau_d_list[i];
    auto [bl_dag, subspace_index_dag] = findIndex(block_shape, iota_d_dag_list[i]);
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
  if (bl_indx == -1) { //-1 is for returning the norm
    auto hyb_mat = diagram::hyb_matrix_t(diagram, Delta_tau);
    sign         = diagram.sign();
    hyb_weight   = inclusion_exclusion(diagram, hyb_mat);
    return hyb_weight * sign * norm(u_products);
  } else if (u_products[bl_indx].size() != 0) {
    auto hyb_mat = diagram::hyb_matrix_t(diagram, Delta_tau);
    sign         = diagram.sign();
    hyb_weight   = inclusion_exclusion(diagram, hyb_mat);
    int bl_size  = std::sqrt(u_products[bl_indx].size());
    int i        = subspace_indx / bl_size;
    int j        = subspace_indx % bl_size;
    return u_products[bl_indx](i, j) * hyb_weight * sign;
  } else {
    return 0.0;
  }
}

template <typename T> inline void print_rank(xfac::TensorTrain<T> tt) {
  int len = tt.M.size();
  std::vector<int> rs(len - 1);
  for (auto i = 0u; i < len - 1; i++) rs[i] = tt.M[i].n_slices;
  std::cout << "rank: ";
  print_vector(rs);
  std::cout << std::endl;
}



