#pragma once
#define ARMA_DONT_USE_OPENMP
#include <armadillo>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <Eigen/Eigen>
#include <unsupported/Eigen/MatrixFunctions>
#include "nvtx.hpp"
#include <random>  
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
#include <inchworm/params.hpp>

using namespace inchworm;

using cv_func = std::function<std::vector<double>(const std::vector<double> &, double, double)>;
using jb_func = std::function<double(const std::vector<double> &, double, double)>;

enum debug_t {
  none, //0, no debug
  low,  //1, simulation level debug
  high  //2, simulation level debug + TCI level debug (print both pivot error and integral)
};

inline std::tuple<int, int, int> bl1_to_bl3(int index, const std::vector<int> &block_shape) {
  int running_sum = 0;
  int bl_indx     = 0;
  int i           = 0;
  int j           = 0;
  for (int bl = 0; bl < block_shape.size(); ++bl) {
    int new_sum = running_sum + block_shape[bl]*block_shape[bl];
    if (new_sum > index) {
      bl_indx            = bl;
      int subspace_index = index - running_sum;
      i                  = subspace_index / block_shape[bl];
      j                  = subspace_index % block_shape[bl];
      return std::make_tuple(bl_indx, i, j);
    }
    running_sum = new_sum;
  }
  std::cerr << "bl1_to_bl3: index out of range\n";
  std::exit(EXIT_FAILURE);
  return std::make_tuple(bl_indx, i, j);
}

inline std::pair<int, int> bl1_to_bl2(int index, const std::vector<int> &block_shape) {
  int running_sum = 0;
  int bl_indx     = 0;
  for (int bl = 0; bl < block_shape.size(); ++bl) {
    int new_sum = running_sum + block_shape[bl]*block_shape[bl];
    if (new_sum > index) {
      int subspace_index = index - running_sum;
      return std::make_pair(bl, subspace_index);
    }
    running_sum = new_sum;
  }
  std::cerr << "bl1_to_bl2: index out of range\n";
  std::exit(EXIT_FAILURE);
  return std::make_pair(bl_indx, index);
}

inline int bl2_to_bl1(int bl_indx, int subspace_index, const std::vector<int> &block_shape) {
  int index = 0;
  for (int bl = 0; bl < bl_indx; ++bl) { index += block_shape[bl]*block_shape[bl]; }
  index += subspace_index;
  return index;
}

inline std::tuple<int, int, int> bl2_to_bl3(int bl_indx, int subspace_index, const std::vector<int> &block_shape) {
  int i = subspace_index / block_shape[bl_indx];
  int j = subspace_index % block_shape[bl_indx];
  return std::make_tuple(bl_indx, i, j);
}

inline std::tuple<int, int> bl3_to_bl1(int bl_indx, int i, int j, const std::vector<int> &block_shape) {
  int index = 0;
  for (int bl = 0; bl < bl_indx; ++bl) { index += block_shape[bl]* block_shape[bl]; }
  index += i * block_shape[bl_indx] + j;
  return std::make_tuple(bl_indx, index);
}

inline std::pair<std::vector<double>, std::vector<double>> generate_inchworm_grid(double ti, double tf, long n_linear, int order_Chebyshev) {
  // n_linear points are the inchworm grid, which has n_linear-1 intervals
  // order_Chebyshev is the order of the Chebyshev approximation within each interval, i.e., order_Chebyshev+1 points are used in each interval
  std::vector<double> grid_linear;
  std::vector<double> grid;
  for (int i = 0; i < n_linear; i++) {
    double wr = static_cast<double>(i) / (n_linear - 1);
    grid_linear.push_back(ti * (1 - wr) + tf * wr);
  }
  for (int i = 0; i < n_linear - 1; i++) {
    double a = grid_linear[i];
    double b = grid_linear[i + 1];
    grid.push_back(a);
    for (int j = order_Chebyshev; j >= 0; j--) {
      double x = 0.5 * (a + b) + 0.5 * (b - a) * cos(M_PI * (2 * j + 1) / (2 * (order_Chebyshev + 1)));
      grid.push_back(x);
    }
  }
  grid.push_back(tf);
  return {grid_linear, grid};
}

inline double get_random_number() {
  // Create a random number generator
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(-1.0, 1.0); // Define the distribution between -1 and 1
  // Generate and return a random number
  return dis(gen);
}

inline double get_hash_random_number(const std::vector<double> &vec) {
  std::size_t seed = 0;
  for (double i : vec) {
    // Hash individual double and combine
    seed ^= std::hash<double>{}(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }

  // Map the resulting hash to a double in [-1, 1]
  // static_cast is used to convert size_t to double
  // This mapping is a simple example, and there are many ways to do it.
  return static_cast<double>(seed) / static_cast<double>(std::numeric_limits<std::size_t>::max()) * 2.0 - 1.0;
}

//printing helper
template <typename T> void print_vector(const std::vector<T> &vec) {
  for (auto v : vec) std::cout << v << ' ';
  std::cout << std::endl;
}

template <typename T> inline void print_rank(T tt) {
  int len = tt.M.size();
  std::vector<int> rs(len - 1);
  for (auto i = 0u; i < len - 1; i++) rs[i] = tt.M[i].n_slices;
  std::cout << "rank: ";
  print_vector(rs);
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
    case 50: return quadrature_GK<50>(a, b);
    case 31: return quadrature_GK<31>(a, b);
    case 29: return quadrature_GK<29>(a, b);
    case 27: return quadrature_GK<27>(a, b);
    case 25: return quadrature_GK<25>(a, b);
    case 23: return quadrature_GK<23>(a, b);
    case 21: return quadrature_GK<21>(a, b);
    case 17: return quadrature_GK<17>(a, b);
    case 19: return quadrature_GK<19>(a, b);
    case 13: return quadrature_GK<13>(a, b);
    default: {
      std::cerr << "select_quadrature_GK: value not supported\n";
      std::exit(EXIT_FAILURE);
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
      std::cerr << "get_elements: index out of range\n";
      std::exit(EXIT_FAILURE);
    }
  }
  return result;
}

template <typename T> std::vector<int> get_elements_int(const std::vector<int> &indices, const std::vector<T> &values) {
  std::vector<int> result;
  for (int index : indices) {
    if (index >= 0 && index < values.size()) {
      result.push_back(static_cast<int>(values[index]));
    } else {
      std::cerr << "get_elements: index out of range\n";
      std::exit(EXIT_FAILURE);
    }
  }
  return result;
}

inline std::vector<double> change_variable0(const std::vector<double> &nus, double tau_max, double tau_min = 0.0) {
  std::vector<double> taus(nus.size());
  taus[0] = nus[0] * (tau_max - tau_min) + tau_min;

  for (size_t it = 1; it < taus.size(); ++it) { taus[it] = taus[it - 1] + nus[it] * (tau_max - taus[it - 1]); }

  return taus;
}

inline double jacobian0(const std::vector<double> &taus, double tau_max, double tau_min = 0.0) {
  double prod = tau_max - tau_min;
  for (size_t j = 1; j < taus.size(); ++j) { prod *= (tau_max - taus[j - 1]); }
  return std::abs(prod);
}

inline std::vector<double> change_variable1(const std::vector<double> &nus, double tau_max, double tau_min = 0.0) {
  std::vector<double> taus(nus.size());
  taus[nus.size() - 1] = tau_min + (tau_max - tau_min) * std::pow((nus[nus.size() - 1]), (1.0 / nus.size()));
  for (int i = nus.size() - 2; i >= 0; --i) { taus[i] = tau_min + (taus[i + 1] - tau_min) * std::pow((nus[i]), (1.0 / (i + 1))); }
  return taus;
}

inline double jacobian1(const std::vector<double> &taus, double tau_max, double tau_min = 0.0) {
  return std::pow((tau_max - tau_min), (taus.size())) / factorial(taus.size());
}

inline std::vector<double> change_variable2(const std::vector<double> &nus, double tau_max, double tau_min = 0.0) {
  std::vector<double> taus(nus.size());
  taus[nus.size() - 1] = tau_max - (tau_max - tau_min) * std::pow((nus[nus.size() - 1]), (1.0 / nus.size()));
  for (int i = nus.size() - 2; i >= 0; --i) { taus[i] = tau_max - (tau_max - taus[i + 1]) * std::pow((nus[i]), (1.0 / (i + 1))); }
  return taus;
}

inline double jacobian2(const std::vector<double> &taus, double tau_max, double tau_min = 0.0) {
  return std::pow((tau_max - tau_min), (taus.size())) / factorial(taus.size());
}

typedef Eigen::Matrix<double, -1, -1, Eigen::ColMajor> DColMatrix;
typedef Eigen::Matrix<double, -1, 1, Eigen::ColMajor> DColVector;
typedef Eigen::Matrix<double, -1, -1, Eigen::RowMajor> DMatrix;

inline DMatrix get_P(unsigned k, double max_val, double min_val) {
  DMatrix P = DColMatrix::Constant(k, k + 1, min_val);

  for (auto i = 0u; i < k; i++) P(i, 0) = max_val;

  for (auto j = 2u; j < k + 1; j++) {
    for (auto i = j - 1; i < k; i++) { P(i, j) = max_val; }
  }
  return P;
}

inline std::vector<double> change_variable3(const std::vector<double> &vs, double max_val_target, double min_val_target = 0.0) {
  double max_val_source = 1.0;
  double min_val_source = 0.0;
  double scale_source   = max_val_source - min_val_source;
  auto P                = get_P(vs.size(), max_val_target, min_val_target);

  auto t = P.col(0);
  for (auto i = 1u; i <= vs.size(); i++) {
    double x = (vs[i - 1] - min_val_source) / scale_source;
    t        = std::pow((x), (1.0 / i)) * t + (1 - std::pow((x), (1.0 / i))) * P.col(i);
  }

  std::vector<double> result(t.size());
  for (int i = 0; i < t.size(); ++i) { result[i] = t(i); }
  return result;
}

inline double jacobian3(const std::vector<double> &vs, double max_val_target, double min_val_target = 0.0) {
  double max_val_source = 1.0;
  double min_val_source = 0.0;
  double scale_source   = max_val_source - min_val_source;
  double scale_target   = max_val_target - min_val_target;

  int k = vs.size();

  return std::pow((scale_target / scale_source), (k)) / std::tgamma(k + 1);
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
                                 int bl_indx, int subspace_indx, bool use_bare_propagator, std::vector<int> const &gf_index,
                                 gf_struct_t const &gf_struct, double energy_shift = 0.0) {
  
  NVTX_RANGE("evaluate_u_tau_max", 0);
  auto config = config_t(frame_zeroth_order, cp.gf_struct, {0.0, tau_split});
  for (auto i : range(tau_d_list.size())) {
    NVTX_RANGE("evaluate construction", 3);
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

  auto diagram = diagram::time_diagram_t{config, {tau_split}};
  frame_t u_products;
  if (gf_index.size() == 0) { // partiion function or propagator
    if (!use_bare_propagator) {
      NVTX_RANGE("impurity product", 5);
      u_products = make_frame(impurity_product(ad_imp, diagram, tau_max, tau_split, energy_shift, &u_interpolator)
                              * impurity_product(ad_imp, diagram, tau_split, 0, energy_shift, &u_interpolator));
    } else {
      u_products = make_frame(impurity_product(ad_imp, diagram, tau_max, tau_split, energy_shift));
    }
    int sign          = 0;
    double hyb_weight = 0.0;
    if (use_bare_propagator) {
      if (bl_indx == -1) { //-1 is for returning the trace
        if (has_zero_trace(ad_imp, diagram)) { return 0.0; }
        auto hyb_mat = diagram::hyb_matrix_t(diagram, Delta_tau);
        sign         = diagram.sign();
        hyb_weight   = hyb_mat.det();
        return hyb_weight * sign * trace(u_products);
      } else if (u_products[bl_indx].size() != 0) {
        NVTX_RANGE("hyb det", 6);
        auto hyb_mat = diagram::hyb_matrix_t(diagram, Delta_tau);
        sign         = diagram.sign();
        hyb_weight   = hyb_mat.det();
        int bl_size  = std::sqrt(u_products[bl_indx].size());
        EXPECTS(bl_size == ad_imp.get_subspace_dim(bl_indx))
        int i = subspace_indx / bl_size;
        int j = subspace_indx % bl_size;
        return u_products[bl_indx](i, j) * hyb_weight * sign;
      } else {
        return 0.0;
      }
    } // end of if (use_bare_propagator)
    else if (u_products[bl_indx].size() != 0) {
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
  } else if (gf_index.size() == 2) { // greens function
    // if (has_zero_trace(ad_imp, diagram)) { return 0.0; }
    auto l                          = impurity_product(ad_imp, diagram, tau_max, tau_split, energy_shift, &u_interpolator);
    auto r                          = impurity_product(ad_imp, diagram, tau_split, 0, energy_shift, &u_interpolator);
    auto [i, bl_d]                  = findIndex(block_shape, gf_index[0]);
    auto [j, bl_dag]                = findIndex(block_shape, gf_index[1]);
    auto [bl_name_d, bl_size_d]     = gf_struct[bl_d];
    auto [bl_name_dag, bl_size_dag] = gf_struct[bl_dag];
    auto l_x_di                     = l * get_op_block_matrix(ad_imp, bl_name_d, i, false);
    auto r_x_djdag                  = r * get_op_block_matrix(ad_imp, bl_name_dag, j, true);
    auto prod                       = make_frame(l_x_di * r_x_djdag);
    auto hyb_mat                    = diagram::hyb_matrix_t(diagram, Delta_tau);
    int sign                        = diagram.sign();
    double hyb_weight               = inclusion_exclusion(diagram, hyb_mat);
    // std::cout << "hyb_weight: " << hyb_weight << std::endl;
    // std::cout << "sign: " << sign << std::endl;
    // std::cout << "trace: " << trace(prod) << std::endl;
    // std::cout << "trace without ds: " << trace(make_frame(l * r)) << std::endl;
    return -1 * trace(prod) * hyb_weight * sign;
  }
  return 0.0;
}

// tci helper

inline void generateCombinationsHelper(const std::vector<int> &n_index_list, std::vector<int> &current_combination, int index,
                                       std::vector<std::vector<int>> &result) {
  if (index == n_index_list.size()) {
    result.push_back(current_combination);
    return;
  }

  for (int i = 0; i < n_index_list[index]; ++i) {
    current_combination[index] = i;
    generateCombinationsHelper(n_index_list, current_combination, index + 1, result);
  }
}

inline std::vector<std::vector<int>> generateCombinations(const std::vector<int> &n_index_list) {
  std::vector<std::vector<int>> result;
  std::vector<int> current_combination(n_index_list.size(), 0);
  generateCombinationsHelper(n_index_list, current_combination, 0, result);
  return result;
}

template <class T, class Index>
double get_integral_ctt_GK(xfac::CTensorTrain<T, Index> const &ctt, std::vector<std::pair<double, double>> const &bounds, int max_depth = 3,
                           double tol = 1e-14) {
  auto M = ctt.M;
  arma::Mat<T> prod(1, 1, arma::fill::eye);
  for (int i = 0; i < M.size(); i++) {
    auto M_i    = M[i];
    auto M_i_0  = M_i(1.0);
    auto [a, b] = bounds[i];
    arma::mat M_i_int(M_i_0.n_rows, M_i_0.n_cols);
// #pragma omp parallel for collapse(2)
    for (size_t ii = 0; ii < M_i_0.n_rows; ++ii) {
      for (size_t jj = 0; jj < M_i_0.n_cols; ++jj) {
        auto func = [&](double x) {
          auto M_i_x = M_i(x);
          return M_i_x.at(ii, jj);
        };
        double error;
        auto integral      = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(func, a, b, max_depth, tol, &error);
        M_i_int.at(ii, jj) = integral;
      }
    }
    prod = prod * M_i_int;
  }
  return prod.eval()(0, 0);
}

template <class T, class Index>
double get_integral_ctt_tanh_sinh(xfac::CTensorTrain<T, Index> const &ctt, std::vector<std::pair<double, double>> const &bounds, double tol) {
  auto M = ctt.M;
  arma::Mat<T> prod(1, 1, arma::fill::eye);
  for (int i = 0; i < M.size(); i++) {
    auto M_i    = M[i];
    auto M_i_0  = M_i(1.0);
    auto [a, b] = bounds[i];
    arma::mat M_i_int(M_i_0.n_rows, M_i_0.n_cols);
// #pragma omp parallel for collapse(2)
    for (size_t ii = 0; ii < M_i_0.n_rows; ++ii) {
      for (size_t jj = 0; jj < M_i_0.n_cols; ++jj) {
        auto func = [&](double x) {
          auto M_i_x = M_i(x);
          return M_i_x.at(ii, jj);
        };
        boost::math::quadrature::tanh_sinh<double> integrator;
        double termination = std::sqrt(tol);
        auto integral      = integrator.integrate(func, a, b, termination);
        M_i_int.at(ii, jj) = integral;
      }
    }
    prod = prod * M_i_int;
  }
  return prod.eval()(0, 0);
}

inline std::pair<std::vector<double>, std::vector<double>> tanh_sinh_quadrature(double a, double b, double h, int k) {

  std::vector<double> abscissas(2 * k + 1);
  std::vector<double> weights(2 * k + 1);

  // Calculate abscissas
  for (int j = 0; j <= 2 * k; ++j) {
    int i        = j - k;
    double t     = i * h;
    double z     = M_PI_2 * sinh(t); // M_PI_2 for pi/2
    abscissas[j] = (b + a) / 2.0 + (b - a) / 2.0 * tanh(z);
  }

  // Calculate weights
  double scaling_factor = (b - a) / 2.0 * M_PI_2 * h;
  for (int j = 0; j <= 2 * k; ++j) {
    int i      = j - k;
    double t   = i * h;
    double z   = M_PI_2 * sinh(t);
    weights[j] = scaling_factor * cosh(t) / (cosh(z) * cosh(z));
  }

  return {abscissas, weights};
}

template <class T> std::vector<T> partial_integral_tt(xfac::TensorTrain<T> const &tt, std::vector<std::vector<double>> const &weight, size_t n_skip) {
  std::vector<T> res{};
  if (n_skip > weight.size()) {
    std::cerr << "n_skip is larger than the number of physical indices." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  auto tt_sum = xfac::TT_sum(tt, weight);
  // index 0 to n_skip-1 are the physical indices that is not integrated over.
  if (n_skip == 0) { return {arma::dot(tt_sum.L[1], tt_sum.R[0])}; }
  auto right_integral = tt_sum.R[n_skip - 1];
  std::vector<int> n_index_list;
  for (size_t i = 0; i < n_skip; i++) { n_index_list.push_back(weight[i].size()); }
  // std::cout<<"n_index_list: ";
  // print_vector(n_index_list);
  auto index_combinations = generateCombinations(n_index_list);
  // std::cout<<"index_combinations: ";
  // for (auto index_combination : index_combinations) {
  //   print_vector(index_combination);
  // }
  for (auto index_combination : index_combinations) {
    auto partial_integral = right_integral;
    for (size_t i = n_skip - 1; i > 0; i--) {
      arma::Mat<T> current_mat = tt.M[i].col(index_combination[i]);
      partial_integral         = current_mat * partial_integral;
    }
    arma::Mat<T> current_mat = tt.M[0].col(index_combination[0]);
    partial_integral         = current_mat * partial_integral;
    res.push_back(partial_integral.eval()(0));
  }
  return res;
}

// check if there are identical elements in a vector
inline bool is_duplicated(const std::vector<double> &taus) {
  std::vector<double> taus_copy = taus;
  std::sort(taus_copy.begin(), taus_copy.end());
  auto last = std::unique(taus_copy.begin(), taus_copy.end());
  return last != taus_copy.end();
}

inline bool if_contains(const std::vector<double> &taus, double tau) {
  return std::find(taus.begin(), taus.end(), tau) != taus.end();
}

// template <typename T_input, typename T_output>
// T_output tci_error_integral(std::function<T_output(std::vector<T_input>)> f, xfac::TensorTrain<T_output> tt, std::vector<std::vector<T_input>> input,
//                             size_t numEval = 1e3) {

//   T_output e = 0; // Error
//   T_output m = 0; // Magnitude
//   std::random_device rd;
//   std::mt19937 mt(rd());
//   std::vector<int> idxs(input.size(), 0);
//   std::vector<T_input> inputs(input.size(), 0);
//   for (size_t sample = 0; sample < numEval; sample++) {
//     for (auto i = 0u; i < idxs.size(); i++) {
//       int local_dim = input[i].size();
//       idxs[i]       = mt() % (local_dim);
//       inputs[i]     = input[i][idxs[i]];
//     }
//     T_output tt_res    = tt.eval(idxs);
//     T_output current_f = f(inputs);
//     m += std::abs(current_f);
//     e += std::abs(tt_res - current_f);
//   }
//   return e / m;
// }

template <typename T_input, typename T_output>
T_output tci_error_integral(std::function<T_output(std::vector<T_input>)> f, xfac::TensorTrain<T_output> tt, std::vector<std::vector<T_input>> input,
                            size_t numEval = 1e3) {

  T_output e = 0; // Error
  T_output m = 0; // Magnitude

// OpenMP for parallelization:
// #pragma omp parallel reduction(+ : e, m)
  {
    std::random_device rd;
    std::mt19937 mt(rd()); // Each thread should have its own random generator
    std::vector<int> idxs(input.size(), 0);
    std::vector<T_input> inputs(input.size(), 0);

// #pragma omp for // Distribute iterations of the outer loop
    for (size_t sample = 0; sample < numEval; sample++) {
      for (auto i = 0u; i < idxs.size(); i++) {
        int local_dim = input[i].size();
        idxs[i]       = mt() % (local_dim);
        inputs[i]     = input[i][idxs[i]];
      }
      T_output tt_res    = tt.eval(idxs);
      T_output current_f = f(inputs);

      // Accumulate e and m with local thread-safe updates
      e += std::abs(tt_res - current_f);
      m += std::abs(current_f);
    }
  } // End of OpenMP parallel region

  return e / m;
}

// template <typename T_input, typename T_output>
// T_output tci_error_integrand(std::function<T_output(std::vector<T_input>)> f, xfac::TensorTrain<T_output> tt, std::vector<std::vector<T_input>> input,
//                              size_t numEval = 1e3) {
//   std::random_device rd;
//   std::mt19937 mt(rd());
//   std::vector<int> idxs(input.size(), 0);
//   std::vector<T_input> inputs(input.size(), 0);
//   double max_error = 0;
//   for (size_t sample = 0; sample < numEval; sample++) {
//     for (auto i = 0u; i < idxs.size(); i++) {
//       int local_dim = input[i].size();
//       idxs[i]       = mt() % (local_dim);
//       inputs[i]     = input[i][idxs[i]];
//     }
//     T_output tt_res    = tt.eval(idxs);
//     T_output current_f = f(inputs);
//     max_error = std::max(max_error, std::abs(tt_res - current_f));
//   }
//   return max_error;
// }

template <typename T_input, typename T_output>
T_output tci_error_integrand(std::function<T_output(std::vector<T_input>)> f, xfac::TensorTrain<T_output> tt, std::vector<std::vector<T_input>> input,
                             size_t numEval = 1e3) {

  double max_error = 0;

// OpenMP for parallelization:
// #pragma omp parallel
  {
    std::random_device rd;
    std::mt19937 mt(rd()); // Each thread should have its own random generator
    std::vector<int> idxs(input.size(), 0);
    std::vector<T_input> inputs(input.size(), 0);
    double thread_max_error = 0; // Local max_error for each thread

// #pragma omp for // Distribute iterations of the outer loop
    for (size_t sample = 0; sample < numEval; sample++) {
      for (auto i = 0u; i < idxs.size(); i++) {
        int local_dim = input[i].size();
        idxs[i]       = mt() % (local_dim);
        inputs[i]     = input[i][idxs[i]];
      }
      T_output tt_res    = tt.eval(idxs);
      T_output current_f = f(inputs);
      thread_max_error   = std::max(thread_max_error, std::abs(tt_res - current_f));
    }

// Critical section to update global max_error safely
// #pragma omp critical
    max_error = std::max(max_error, thread_max_error);
  } // End of OpenMP parallel region

  return max_error;
}

inline double quadrature_1D_ts(std::function<double(double)> const &func, double a, double b) {
  int n_grid_tanh_sinh  = 15;
  double h_tanh_sinh    = 4.0 / n_grid_tanh_sinh;
  auto [xi_old, wi_old] = tanh_sinh_quadrature(a, b, h_tanh_sinh, n_grid_tanh_sinh);
  // find all 0 and 1 in xi and remove the corresponding xi and wi
  std::vector<double> xi, wi;
  double tol = 1e-14;
  for (int i = 0; i < xi_old.size(); i++) {
    if (abs(xi_old[i]) > tol && abs(xi_old[i] - 1) > tol) {
      xi.push_back(xi_old[i]);
      wi.push_back(wi_old[i]);
    }
  }
  double integral_val = 0.0;
  for (int i = 0; i < xi.size(); i++) { integral_val += wi[i] * func(xi[i]); }
  return integral_val;
}

template <class T, class Index>
double get_integral_ctt_tanh_sinh_depth0(xfac::CTensorTrain<T, Index> const &ctt, std::vector<std::pair<double, double>> const &bounds) {
  auto M = ctt.M;
  arma::Mat<T> prod(1, 1, arma::fill::eye);
  for (int i = 0; i < M.size(); i++) {
    auto M_i    = M[i];
    auto M_i_0  = M_i(1.0);
    auto [a, b] = bounds[i];
    arma::mat M_i_int(M_i_0.n_rows, M_i_0.n_cols);
// #pragma omp parallel for collapse(2)
    for (size_t ii = 0; ii < M_i_0.n_rows; ++ii) {
      for (size_t jj = 0; jj < M_i_0.n_cols; ++jj) {
        auto func = [&](double x) {
          auto M_i_x = M_i(x);
          return M_i_x.at(ii, jj);
        };
        auto integral      = quadrature_1D_ts(func, a, b);
        M_i_int.at(ii, jj) = integral;
      }
    }
    prod = prod * M_i_int;
  }
  return prod.eval()(0, 0);
}

template <typename T_output, typename T_input>
T_output calculate_sum_recursion(std::function<T_output(std::vector<T_input>)> integrand, const std::vector<std::vector<T_input>> &choices,
                                 const std::vector<std::vector<double>> &weights, double const_jacobian) {
  T_output sum = 0.0;
  std::vector<int> indices(choices.size(), 0);
  std::vector<T_input> currentChoice(choices.size());
  std::function<void(size_t)> iterateChoices = [&](size_t index) {
    if (index == choices.size()) {
      double product = 1.0;
      for (size_t i = 0; i < choices.size(); ++i) {
        currentChoice[i] = choices[i][indices[i]];
        product *= weights[i][indices[i]];
      }
      sum += integrand(currentChoice) * product;
      return;
    }

    for (size_t i = 0; i < choices[index].size(); ++i) {
      indices[index] = static_cast<int>(i);
      iterateChoices(index + 1);
    }
  };

  {
    iterateChoices(0);
  }

  return sum * const_jacobian;
}

template <typename T_input>
std::pair<std::vector<std::vector<T_input>>, std::vector<double>> generate_combination_at_fixed_order(int order, std::vector<T_input> const &choices,
                                                                                                      std::vector<double> const &weights) {
  std::vector<std::vector<T_input>> res;
  std::vector<double> res_weights;
  std::vector<int> indices(2 * order, 0);
  std::vector<T_input> currentChoice(2 * order);
  std::function<void(size_t)> iterateChoices = [&](size_t index) {
    if (index == 2 * order) {
      double product = 1.0;
      for (size_t i = 0; i < 2 * order; ++i) {
        currentChoice[i] = choices[indices[i]];
        product *= weights[indices[i]];
      }
      res.push_back(currentChoice);
      res_weights.push_back(product);
      return;
    }

    for (size_t i = 0; i < choices.size(); ++i) {
      indices[index] = static_cast<int>(i);
      iterateChoices(index + 1);
    }
  };

  iterateChoices(0);
  return {res, res_weights};
}

template <typename T_output, typename T_input>
std::vector<T_output> calculate_sum(std::function<T_output(std::vector<T_input>)> integrand, const std::vector<std::vector<T_input>> &inputs,
                                    const std::vector<double> &weights, double const_jacobian) {
  T_output sum = 0.0;
// #pragma omp parallel for reduction(+ : sum)
  for (size_t i = 0; i < inputs.size(); i++) { sum += integrand(inputs[i]) * weights[i]; }
  return {sum * const_jacobian};
}

template <typename T_output, typename T_input>
std::vector<T_output> do_TCI(std::function<T_output(std::vector<T_input>)> integrand, std::vector<std::vector<T_input>> const &input,
                             std::vector<std::vector<double>> const &weight, std::vector<int> const &pivot1, long &count, int sweep_bound,
                             int bond_dim, double reltol, bool fullPiv, int tci_prrlu, int error_type, size_t error_eval, double convergence_bound,
                             int convergence_iter, debug_t debug, std::vector<std::vector<int>> const &init_global_pivots, double const_jacobian,
                             int unsummed_tci, int unsummed_tot_size, bool adaptive_error, double decay_rate, double *auxi_height = nullptr) {
  NVTX_RANGE("do_TCI", 1);
  double last_error{0};
  double current_error{0};
  std::vector<T_output> integral{};
  std::vector<T_output> previous_integral{};
  std::vector<std::vector<std::vector<int>>> pivots{};
  std::vector<int> valid_init_global_pivot{};
  if (unsummed_tci !=0 ) {
    integral          = std::vector<T_output>(unsummed_tot_size, 0);
    previous_integral = std::vector<T_output>(unsummed_tot_size, 0);
  } else {
    integral          = std::vector<T_output>(1, 0);
    previous_integral = std::vector<T_output>(1, 0);
  }
  double first_error{0};
  if (tci_prrlu == 1 || tci_prrlu == 2) {
    if (debug > 1) {
      std::cout << "TCI 2" << std::endl;
      std::cout << "iteration nEval error integral\n";
    }
    int bond_dim_init = bond_dim;
    if (tci_prrlu == 2) { bond_dim_init = 1; }
    for (int i = 1; i <= sweep_bound; i++) {
      auto ci =
         xfac::CTensorCI2<T_output, T_input>(integrand, input, {.bondDim = bond_dim_init, .reltol = reltol, .pivot1 = pivot1, .fullPiv = fullPiv});
      if (tci_prrlu == 2) {
        ci.param.bondDim = bond_dim_init + i;
        if (ci.param.bondDim > bond_dim) { ci.param.bondDim = bond_dim; }
      }
      {
        NVTX_RANGE("ci: paste pivots", 6);
      if (i > 1) {
        for (auto b = 0u; b < ci.len() - 1; b++) ci.addPivotsAt(pivots[b], b);
        // ci.iterate(1, 0);
        ci.makeCanonical();
      }
      }
      if (init_global_pivots.size() != 0 && i == 1) { ci.addPivotsAllBonds(init_global_pivots); }
      {
        NVTX_RANGE("ci: iterate", 7);
      ci.iterate();
      }
      if (i == 1) { *auxi_height = ci.pivotError[ci.pivotError.size() - 1]; }
      // ci.makeCanonical();
      auto tt = ci.tt;
      {
        NVTX_RANGE("ci: obtain integral", 9);
      if (unsummed_tci != 0) {
        integral = partial_integral_tt(tt, weight, unsummed_tci);
      } else {
        integral = partial_integral_tt(tt, weight, 0);
      }
      }
      for (auto i = 0; i < integral.size(); i++) { integral[i] *= const_jacobian; }
      // integral = ci.tt.sum(weight);
      {
        NVTX_RANGE("ci: error calculation", 8);
      if (error_type == 0) {
        current_error = ci.pivotError[ci.pivotError.size() - 1];
      } else if (error_type == 1) {
        current_error = ci.trueError(error_eval);
      } else if (error_type == 2) {
        current_error = tci_error_integral(integrand, ci.tt, input, error_eval);
      } else if (error_type == 3) {
        current_error = tci_error_integrand(integrand, ci.tt, input, error_eval);
      } else {
        std::cerr << "error_type not supported" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      }
      if (i == 1) { first_error = current_error; }
      if (debug > 1) { std::cout << i << " " << count << " " << current_error << " " << integral[0] << std::endl; }
      // if (std::abs(current_error / first_error) < convergence_bound && i > convergence_iter) { break; }
      T_output diff = 0;
      for (auto i = 0u; i < integral.size(); i++) { diff += std::abs(previous_integral[i] - integral[i]); }
      if (adaptive_error) {
        std::cout << "auxi_height: " << *auxi_height << std::endl;
        // *auxi_height = *auxi_height * decay_rate;
        // if (last_error < current_error) { *auxi_height = current_error * decay_rate; }
        *auxi_height = ci.pivotError[ci.pivotError.size() - 1] * decay_rate;
      }
      if (diff < convergence_bound && i > convergence_iter) { break; }
      last_error        = current_error;
      previous_integral = integral;
      if (debug > 1) { print_rank(ci.tt); }
      pivots.clear();
      {
        NVTX_RANGE("ci: copy pivots", 8);
      for (auto b = 0u; b < ci.len() - 1; b++) {
        auto pivots_at_b = ci.getPivotsAt(b);
        pivots.push_back(pivots_at_b);
        if (valid_init_global_pivot.size() == 0) {
          // temporaliy change auxi_height to 0
          double current_auxi_height = *auxi_height;
          *auxi_height               = 0;
          for (auto pivot : pivots_at_b) {
            std::vector<T_input> inputs;
            for (auto i = 0u; i < input.size(); i++) { inputs.push_back(input[i][pivot[i]]); }
            if (integrand(inputs) != 0) { valid_init_global_pivot = pivot; }
          }
          *auxi_height = current_auxi_height;
        }
      }
    }
      if (adaptive_error && *auxi_height > convergence_bound && i == sweep_bound) {
        std::cout << "Warning: the last auxiliary height is larger than the convergence bound." << std::endl;
        std::cout << "auxi_height: " << *auxi_height << std::endl;
      }
    }
    if (valid_init_global_pivot.size() == 0) {
      if (unsummed_tci != 0) {
        integral = std::vector<T_output>(unsummed_tot_size, 0);
      } else {
        integral = std::vector<T_output>(1, 0);
      }
      std::cout << "Warning: no valid initial global pivot found." << std::endl;
      return integral;
    }
    *auxi_height = 0;
    auto ci      = xfac::CTensorCI2<T_output, T_input>(
       integrand, input, {.bondDim = bond_dim_init + sweep_bound, .reltol = reltol, .pivot1 = valid_init_global_pivot, .fullPiv = fullPiv});
    for (auto b = 0u; b < ci.len() - 1; b++) ci.addPivotsAt(pivots[b], b);
    // ci.iterate(1, 0);
    ci.makeCanonical();
    auto tt = ci.tt;
    if (unsummed_tci != 0) {
      integral = partial_integral_tt(tt, weight, unsummed_tci);
    } else {
      integral = partial_integral_tt(tt, weight, 0);
    }
    for (auto i = 0; i < integral.size(); i++) { integral[i] *= const_jacobian; }

    // ci.makeCanonical();
    // auto ctt = ci.get_CTensorTrain();
    // std::cout << "ctt obtained" << std::endl;
    // std::vector<std::pair<double, double>> bounds;
    // for (auto i = 0u; i < input.size(); i++) { bounds.push_back({0, 1}); }
    // double integral_val = get_integral_ctt_GK(ctt, bounds,2, 1e-14);
    // // double integral_val = get_integral_ctt_tanh_sinh_depth0(ctt, bounds);
    // integral            = integral_val;
    // std::cout << "integral_adaptive: " << integral << std::endl;
  } else {
    std::cerr << "tci_prrlu not supported" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (debug > 1) { std::cout << std::endl; }
  return integral;
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

inline std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
obtain_taus(const std::vector<double> &vs, int n_left, double tau_split, double tau_max, cv_func change_variable) {
  if (n_left == 0 && tau_split == 0.0) {
    std::vector<double> taus(vs.size());
    taus = change_variable(vs, tau_max, 0.0);
    return std::make_tuple(taus, taus, taus);
  }
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
  double x_val = (x[0] + x[1]) / 2.0;
  return std::sin(x_val * M_PI);
}

template <typename T> inline double sin_func_all(std::vector<T> const &vs, std::vector<T> const &iotas, int n_opts) {
  std::vector<T> iotas_normalized(iotas.size());
  for (int i = 0; i < iotas.size(); ++i) { iotas_normalized[i] = iotas[i] / n_opts; }
  std::vector<T> v_iota_s_normalized(vs.size() + iotas.size());
  std::copy(vs.begin(), vs.end(), v_iota_s_normalized.begin());
  std::copy(iotas_normalized.begin(), iotas_normalized.end(), v_iota_s_normalized.begin() + vs.size());
  double res = 0.0;
  for (int i = 0; i < v_iota_s_normalized.size() - 1; ++i) {
    std::vector<T> pair = {v_iota_s_normalized[i], v_iota_s_normalized[i + 1]};
    res += sin_func_pair(pair);
  }
  res /= (v_iota_s_normalized.size() - 1);
  return res;
}

// template <typename T> inline double sin_func_pair(std::vector<T> const &x) {
//   double x_val = (x[0]+x[1]) / 2.0;
//   return std::sin(2*x_val * M_PI);
// }

// template <typename T> inline double sin_func_all(std::vector<T> const &vs, std::vector<T> const &iotas, int n_opts){
//   std::vector<T> iotas_normalized(iotas.size());
//   for (int i = 0; i < iotas.size(); ++i) { iotas_normalized[i] = (iotas[i]+1) / (n_opts+1); }
//   std::vector<T> v_iota_s_normalized(vs.size()+iotas.size());
//   std::copy(vs.begin(), vs.end(), v_iota_s_normalized.begin());
//   std::copy(iotas_normalized.begin(), iotas_normalized.end(), v_iota_s_normalized.begin()+vs.size());
//   double res = 0.0;
//   for(int i = 0; i < v_iota_s_normalized.size()-1; ++i) {
//     std::vector<T> pair = {v_iota_s_normalized[i], v_iota_s_normalized[i+1]};
//     res += sin_func_pair(pair);
//   }
//   res /= (v_iota_s_normalized.size()-1);
//   return res;
// }

template <typename T> inline double linear_func_all(std::vector<T> const &iotas, int n_opts) {
  double x_val = 0;
  double base  = 1.0 / (n_opts);
  for (int i = 0; i < iotas.size(); ++i) { x_val += iotas[i] * std::pow(base, i + 1); }
  // return (2.0* (x_val + std::pow(base, iotas.size())/2)-1.0)/2.0;
  return std::sin((x_val + std::pow(base, iotas.size()) / 2) * M_PI);
}

struct myOperator {
  int iota{};
  bool is_dagger{};
};

inline std::vector<std::vector<bool>> generate_all_wavefunctions(int N) {
  std::vector<std::vector<bool>> all_vectors;
  int num_vectors = 1 << N; // Calculate 2^N

  for (int i = 0; i < num_vectors; ++i) {
    std::vector<bool> current_vector(N);
    for (int j = 0; j < N; ++j) { current_vector[j] = (i >> j) & 1; }
    all_vectors.push_back(current_vector);
  }

  return all_vectors;
}

inline bool pairCompare(const std::pair<std::vector<int>, std::vector<int>> &a, const std::pair<std::vector<int>, std::vector<int>> &b) {
  // Compare the two pairs
  return a.first == b.first && a.second == b.second;
}

inline void removeDuplicates(std::vector<std::pair<std::vector<int>, std::vector<int>>> &vec) {
  // Sort the vector to bring duplicates together
  std::sort(vec.begin(), vec.end());

  // Use unique() function to remove duplicates
  auto last = std::unique(vec.begin(), vec.end(), pairCompare);

  // Erase the duplicates
  vec.erase(last, vec.end());
}

inline std::vector<std::pair<std::vector<int>, std::vector<int>>> generate_phi_segment(std::vector<double> const &iotas,
                                                                                       std::vector<int> const &gf_block_shape) {
  // check if gf_block_shape is all ones
  std::vector<std::pair<std::vector<int>, std::vector<int>>> phi_pair_list;
  std::vector<int> ones(gf_block_shape.size(), 1);
  if (gf_block_shape != ones) {
    std::cerr << "generate_phi_segment: gf_block_shape is not all ones\n";
    std::exit(EXIT_FAILURE);
  }
  int n_phi = gf_block_shape.size();
  // generate all possible wavefunction configurations

  auto all_wavefunctions = generate_all_wavefunctions(n_phi);
  for (auto is_created : all_wavefunctions) {
    std::vector<myOperator> operators{};
    for (auto iota : iotas) {
      int iota_int = static_cast<int>(iota);
      if (is_created[iota_int]) {
        operators.push_back({iota_int, false});
        is_created[iota_int] = false;
      } else {
        operators.push_back({iota_int, true});
        is_created[iota_int] = true;
      }
    }
    std::vector<int> phi_d_list;
    std::vector<int> phi_d_dag_list;
    //collect the indices for c and c^dagger
    for (size_t i = 0; i < operators.size(); ++i) {
      if (operators[i].is_dagger) {
        phi_d_dag_list.push_back(i);
      } else {
        phi_d_list.push_back(i);
      }
    }
    phi_pair_list.push_back({phi_d_list, phi_d_dag_list});
  }
  removeDuplicates(phi_pair_list);
  return phi_pair_list;
}