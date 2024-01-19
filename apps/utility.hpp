#pragma once
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

enum debug_t {
  none, //0, no debug
  low,  //1, simulation level debug
  high  //2, simulation level debug + TCI level debug
};

//printing helper
template <typename T> void print_vector(const std::vector<T> &vec) {
  for (auto v : vec) std::cout << v << ' ';
  std::cout << std::endl;
}

template <typename T> inline void print_rank(xfac::TensorTrain<T> tt) {
  int len = tt.M.size();
  std::vector<int> rs(len - 1);
  for (auto i = 0u; i < len - 1; i++) rs[i] = tt.M[i].n_slices;
  std::cout << "rank: ";
  print_vector(rs);
  std::cout << std::endl;
}

template <typename T> void print_block_shape(block_gf<imtime, T> const &x_tau) {
  std::cout << "number of taus: " << x_tau[0].mesh().size() << std::endl;
  std::cout << "number of block: " << x_tau.size() << std::endl;
  for (int bl = 0; bl < x_tau.size(); ++bl) { std::cout << "block: " << bl << " shape: " << x_tau[bl].target_shape() << std::endl; }
}

// quadrature helper
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

// generate helper
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
template <typename T> void generate_combinations(const std::vector<T> &A, std::vector<T> &B, int idx, std::vector<std::vector<T>> &all_combinations) {
  if (idx == B.size()) {
    all_combinations.push_back(B);
    return;
  }

  for (int i = 0; i < A.size(); ++i) {
    B[idx] = A[i];
    generate_combinations(A, B, idx + 1, all_combinations);
  }
}

//similar as generate_combinations but for first element of all_combinations, it only takes values from A_half
template <typename T>
void generate_combinations_symmetrized(const std::vector<T> &A, const std::vector<T> &A_half, std::vector<T> &B, int idx,
                                       std::vector<std::vector<T>> &all_combinations) {
  if (idx == B.size()) {
    all_combinations.push_back(B);
    return;
  }
  if (idx == 0) {
    for (int i = 0; i < A_half.size(); ++i) {
      B[idx] = A_half[i];
      generate_combinations(A, B, idx + 1, all_combinations);
    }
  } else {
    for (int i = 0; i < A.size(); ++i) {
      B[idx] = A[i];
      generate_combinations(A, B, idx + 1, all_combinations);
    }
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
//     taus[i] = tau_min + (taus[i+1]- tau_min) * std::pow(nus[i], (1.0 /static_cast<double>(i + 1)));
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

inline double evaluate_u_tau_max(frame_t &frame_zeroth_order, double tau_split, double tau_max, std::vector<std::vector<fop_t>> const &all_d_ops,
                                 std::vector<std::vector<fop_t>> const &all_d_dag_ops, std::vector<int> const &block_shape, constr_params_t const &cp,
                                 hyb_tau_t const &Delta_tau, atom_diag const &ad_imp, interpolator_t<scalar_t> const &u_interpolator,
                                 auto const &tau_d_list, auto const &tau_d_dag_list, auto const &iota_d_list, auto const &iota_d_dag_list,
                                 int bl_indx, int subspace_indx) {
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

template <typename T_output, typename T_input>
T_output do_TCI(std::function<T_output(std::vector<T_input>)> func, std::vector<std::vector<T_input>> &input,
                std::vector<std::vector<double>> &weight, std::vector<int> &pivot1, int sweep_bound, int bond_dim, double integral_error_bound,
                double pivot_error_bound, bool tci_prrlu, debug_t debug, long &count) {
  double current_integral{0};
  double previous_integral{0};
  double last_pivot_error{0};
  if (debug > 1) { std::cout << "iteration nEval LastSweepPivotError integral\n"; }
  if (tci_prrlu) {
    auto ci = xfac::CTensorCI2<T_output, T_input>(func, input, {.bondDim = bond_dim, .reltol = 1e-18, .pivot1 = pivot1, .fullPiv = true});
    std::cout << "bond_dim: " << ci.param.bondDim << std::endl;
    for (int i = 1; i <= sweep_bound; i++) {
      ci.iterate();
      // ci.makeCanonical();
      current_integral = ci.tt.sum(weight);
      last_pivot_error = ci.pivotError[ci.pivotError.size() - 1];
      // last_pivot_error = ci.trueError();
      if (debug > 1) { std::cout << i << " " << count << " " << last_pivot_error << " " << current_integral << std::endl; }
      //  if ( last_pivot_error < pivot_error_bound && i > 1) { break; }
      // if (std::abs(current_integral - previous_integral) < integral_error_bound || last_pivot_error < pivot_error_bound && i > 1) { break; }
      // if (i == 2) {
      //   std::cout << "add pivots" << std::endl;
      //   for (auto b = 0u; b < ci.len() - 1; b++) {
      //     auto pivots = ci.getPivotsAt(b);
      //     auto first_pivots = pivots[0];
      //     for(auto &p : first_pivots) {
      //       p = (p+15)%30;
      //     }
      //     ci.addPivotsAt(pivots, b);
      //   }
      // }
      previous_integral = current_integral;
    }
    if (debug > 1) { print_rank(ci.tt); }
  } else {
    auto ci = xfac::CTensorCI<T_output, T_input>(func, input, {.reltol = 1e-18, .pivot1 = pivot1, .fullPiv = true, .weight = weight});
    for (int i = 1; i <= sweep_bound+1; i++) {
      ci.iterate();
      current_integral = ci.get_TensorTrain().sum(weight);
      last_pivot_error = ci.pivotError[ci.pivotError.size() - 1];
      // last_pivot_error = ci.trueError();
      if (debug > 1) { std::cout << i-1 << " " << count << " " << last_pivot_error << " " << current_integral << std::endl; }
      // if ( last_pivot_error < pivot_error_bound && i > 1) { break; }
      // if (std::abs(current_integral - previous_integral) < integral_error_bound || last_pivot_error < pivot_error_bound && i > 1) { break; }
      previous_integral = current_integral;
    }
    if (debug > 1) {
      print_rank(ci.get_TensorTrain());
    }
  }
  if (debug > 1) { std::cout << std::endl; }
  return current_integral;
}

template <typename T_output, typename T_input>
T_output do_TCI_add_pivots(std::function<T_output(std::vector<T_input>)> func, std::vector<std::vector<T_input>> &input,
                           std::vector<std::vector<double>> &weight, std::vector<int> &pivot1, int sweep_bound, int bond_dim,
                           double integral_error_bound, double pivot_error_bound, bool tci_prrlu, debug_t debug, long &count,
                           std::vector<std::vector<int>> &valid_pivots) {
  double current_integral{0};
  double previous_integral{0};
  double last_pivot_error{0};
  if (debug > 1) { std::cout << "iteration nEval LastSweepPivotError integral\n"; }
  if (tci_prrlu) {
    std::cout << "test" << std::endl;
    auto ci = xfac::CTensorCI2<T_output, T_input>(func, input, {.bondDim = bond_dim, .reltol = 1e-18, .pivot1 = pivot1, .fullPiv = false});
    std::cout << "bond_dim: " << ci.param.bondDim << std::endl;
    ci.myAddPivotsAllBonds(valid_pivots);
    if (debug > 1) { print_rank(ci.tt); }
    // ci.tt.compressCI(ci.param.reltol, ci.param.bondDim);
    ci.makeCanonical();
    if (debug > 1) { print_rank(ci.tt); }
    for (int i = 1; i <= sweep_bound; i++) {
      ci.iterate();
      ci.makeCanonical();
      current_integral = ci.tt.sum(weight);
      last_pivot_error = ci.pivotError[ci.pivotError.size() - 1];
      if (debug > 1) { std::cout << i-1 << " " << count << " " << last_pivot_error << " " << current_integral << std::endl; }
      previous_integral = current_integral;
      if (debug > 1) { print_rank(ci.tt); }
    }
  } else {
    auto ci = xfac::CTensorCI<T_output, T_input>(func, input, {.reltol = 1e-18, .pivot1 = pivot1, .weight = weight});
    for (int i = 1; i <= sweep_bound; i++) {
      ci.iterate();
      current_integral = ci.get_TensorTrain().sum(weight);
      last_pivot_error = ci.pivotError[ci.pivotError.size() - 1];
      if (debug > 1) { std::cout << i-1 << " " << count << " " << last_pivot_error << " " << current_integral << std::endl; }
      // if ( last_pivot_error < pivot_error_bound && i > 1) { break; }
      // if (std::abs(current_integral - previous_integral) < integral_error_bound || last_pivot_error < pivot_error_bound && i > 1) { break; }
      previous_integral = current_integral;
    }
    if (debug > 1) {
      std::cout << "rank:" << std::endl;
      print_rank(ci.get_TensorTrain());
    }
  }
  if (debug > 1) { std::cout << std::endl; }
  return current_integral;
}

template <typename T_output, typename T_input>
T_output do_TCI_reuse_pivots(std::function<T_output(std::vector<T_input>)> func, std::vector<std::vector<T_input>> &input,
                             std::vector<std::vector<double>> &weight, std::vector<int> &pivot1, int sweep_bound, int bond_dim,
                             double integral_error_bound, double pivot_error_bound, bool tci_prrlu, debug_t debug, long &count,
                             std::vector<std::vector<std::vector<int>>> &previous_pivots) {
  double current_integral{0};
  double previous_integral{0};
  double last_pivot_error{0};
  if (debug > 1) { std::cout << "iteration nEval LastSweepPivotError integral\n"; }
  if (tci_prrlu) {
    auto ci = xfac::CTensorCI2<T_output, T_input>(func, input, {.bondDim = bond_dim, .pivot1 = pivot1});
    //only supported in prrlu
    if (!previous_pivots.empty()) {
      if (debug > 1) { std::cout << "reuse pivots" << std::endl; }
      for (auto b = 0u; b < ci.len() - 1; b++) {
        auto pivots = previous_pivots[b];
        // auto first_half_pivots = std::vector(pivots.begin(), pivots.begin() + pivots.size() / 2);
        ci.addPivotsAt(pivots, b);
      }
      // ci.makeCanonical();
      // last_pivot_error = ci.pivotError[ci.pivotError.size() - 1];
      current_integral = ci.tt.sum(weight);
      // if (std::abs(current_integral - previous_integral) > integral_error_bound) {
      //   previous_integral = current_integral;
      //   for (int i = 1; i <= sweep_bound; i++) {
      //     ci.iterate();
      //     ci.makeCanonical();
      //     last_pivot_error = ci.pivotError[ci.pivotError.size() - 1];
      //     current_integral = ci.tt.sum(weight);
      //     if (debug > 1) { std::cout << i << " " << count << " " << last_pivot_error << " " << current_integral << std::endl; }
      //     if (std::abs(current_integral - previous_integral) < integral_error_bound || last_pivot_error < pivot_error_bound) { break; }
      //     previous_integral = current_integral;
      //   }
      // }
    } else { //empty
      for (int i = 1; i <= sweep_bound; i++) {
        ci.iterate();
        ci.makeCanonical();
        current_integral = ci.tt.sum(weight);
        last_pivot_error = ci.pivotError[ci.pivotError.size() - 1];
        if (debug > 1) { std::cout << i << " " << count << " " << last_pivot_error << " " << current_integral << std::endl; }
        if (std::abs(current_integral - previous_integral) < integral_error_bound || last_pivot_error < pivot_error_bound) { break; }
        previous_integral = current_integral;
      }
    }
    if (debug > 1) { print_rank(ci.tt); }
    if (previous_pivots.empty()) {
      previous_pivots.reserve(ci.len() - 1);
      for (auto b = 0u; b < ci.len() - 1; b++) {
        auto pivots = ci.getPivotsAt(b);
        previous_pivots.push_back(pivots);
      }
    }
  } else {
    auto ci = xfac::CTensorCI<T_output, T_input>(func, input, {.pivot1 = pivot1});
    for (int i = 1; i <= sweep_bound; i++) {
      ci.iterate();
      current_integral = ci.get_TensorTrain().sum(weight);
      last_pivot_error = ci.pivotError[ci.pivotError.size() - 1];
      if (debug > 1) { std::cout << i << " " << count << " " << last_pivot_error << " " << current_integral << std::endl; }
      if (std::abs(current_integral - previous_integral) < integral_error_bound || last_pivot_error < pivot_error_bound) { break; }
      previous_integral = current_integral;
    }
    if (debug > 1) {
      std::cout << "rank:" << std::endl;
      print_rank(ci.get_TensorTrain());
    }
  }
  if (debug > 1) { std::cout << std::endl; }
  return current_integral;
}

inline void print_pivot1(std::vector<int> const &iota_d_list, std::vector<int> const &iota_d_dag_list, std::vector<double> const &tau_d_list,
                         std::vector<double> const &tau_d_dag_list, double pivot_value) {
  std::cout << "iota_d_list: ";
  print_vector(iota_d_list);
  std::cout << "iota_d_dag_list: ";
  print_vector(iota_d_dag_list);
  std::cout << "tau_d_list: ";
  print_vector(tau_d_list);
  std::cout << "tau_d_dag_list: ";
  print_vector(tau_d_dag_list);
  std::cout << "get_u_tau_max_element(pivot1): " << pivot_value << std::endl;
  std::cout << std::endl;
}

inline std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> obtain_taus(const std::vector<double> &vs, int n_left,
                                                                                             double tau_split, double tau_max) {
  std::vector<double> vs_left(vs.begin(), vs.begin() + n_left);
  std::vector<double> vs_right(vs.begin() + n_left, vs.end());
  std::vector<double> taus_left  = change_variable(vs_left, tau_split, 0.0);
  std::vector<double> taus_right = change_variable(vs_right, tau_max, tau_split);
  auto taus(taus_left);
  taus.insert(taus.end(), taus_right.begin(), taus_right.end());

  return std::make_tuple(taus_left, taus_right, taus);
}

template <typename T1, typename T2> void sort_B_according_A(std::vector<T1> &A, std::vector<T2> &B, double reltol = 1e-20) {
  std::vector<size_t> indices(A.size());
  std::iota(indices.begin(), indices.end(), 0); // Fill with 0, 1, 2, ...

  // Sort indices based on corresponding values in A
  std::sort(indices.begin(), indices.end(), [&](size_t i1, size_t i2) { return std::abs(A[i1]) > std::abs(A[i2]); });

  // Create a sorted copy of B and A
  std::vector<T1> sorted_A(A.size());
  std::vector<T2> sorted_B(B.size());
  for (size_t i = 0; i < indices.size(); ++i) {
    sorted_A[i] = A[indices[i]];
    sorted_B[i] = B[indices[i]];
  }
  // truncation according to reltol
  T1 A0                   = sorted_A[0];
  size_t truncation_index = 0;
  for (size_t i = 0; i < sorted_A.size(); ++i) {
    truncation_index = i;
    if (std::abs(sorted_A[i]) < std::abs(A0) * reltol) { break; }
  }
  truncation_index++;
  //print A and sorted A
  // std::cout<< "A:  ";
  // print_vector(A);
  // std::cout<< "sorted_A:  ";
  // print_vector(sorted_A);

  B = std::vector<T2>(sorted_B.begin(), sorted_B.begin() + truncation_index);
  A = std::vector<T1>(sorted_A.begin(), sorted_A.begin() + truncation_index);
}

template <typename T> inline double sin_func(std::vector<T> const &x, int n_opt) {
  double x_val = 0;
  // std::cout << "n_opt: " << n_opt << std::endl;
  // std::cout << "x.size(): " << x.size() << std::endl;
  // compress x between 0 and  1
  double base = 1.0 / (n_opt);
  // std::cout << "base: " << base << std::endl;
  for (int i = 0; i < x.size(); ++i) { x_val += x[i] * std::pow(base, i + 1); }
  // std::cout << "x_val: " << x_val << std::endl;
  return std::sin(x_val * M_PI);
}

template <typename T> inline double sin_func_pair(std::vector<T> const &x) {
  double x_val = (x[0]+x[1]) / 2.0;
  return std::sin(x_val * M_PI);
}

template <typename T> inline double sin_func_all(std::vector<T> const &vs, std::vector<T> const &iotas, int n_opts){
  std::vector<T> iotas_normalized(iotas.size());
  for (int i = 0; i < iotas.size(); ++i) { iotas_normalized[i] = iotas[i] / n_opts; }
  std::vector<T> v_iota_s_normalized(vs.size()+iotas.size());
  std::copy(vs.begin(), vs.end(), v_iota_s_normalized.begin());
  std::copy(iotas_normalized.begin(), iotas_normalized.end(), v_iota_s_normalized.begin()+vs.size());
  double res = 0.0;
  for(int i = 0; i < v_iota_s_normalized.size()-1; ++i) {
    std::vector<T> pair = {v_iota_s_normalized[i], v_iota_s_normalized[i+1]};
    res += sin_func_pair(pair);
  }
  res /= (v_iota_s_normalized.size()-1);
  return res;
}
