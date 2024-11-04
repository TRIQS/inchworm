#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <string>
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
#include "hubbard.hpp"
#include "utility.hpp"

struct global_params_t {
  std::string target{};
  std::string integrand{};
  std::string integral_variable{};
  std::string tci_shape{};
  std::string ergodicity{};
  int model_type{}; // 0 for discrete bath, 1 for continuous bath, 2 for bethe lattice
  bool do_segment = false;
  bool do_cache   = true;
  std::string output_prefix{};
  int unsummed_tci                = 0; // 0: all indices are summed; 1: the first index is not summed; 2: the first two indices are not summed ... 
  double energy_shift              = 0.0;
  double Z_energy_shift_correction = 1.0;
  bool do_adaptive_nGK           = false;
  bool do_global_pivot         = false;
  int map_type = 0;
};

struct model_params_t {
  int n_site{};
  int n_bath{};
  int n_spin{};
  double rescale{};
  double U{};
  double mu{};
  double t{};
  mat_t theta{};
  vec_t epsilon{};
  // the following are derived
  hyb_tau_t Delta_tau{};
  atom_diag ad_imp{};
  std::vector<int> gf_block_shape{};
  std::vector<std::vector<fop_t>> all_d_ops{};
  std::vector<std::vector<fop_t>> all_d_dag_ops{};
  long n_bl{};
  fundamental_operator_set fops{};
  int n_phi{};
};

struct tci_params_t {
  int n_GK{};
  int tci_prrlu{};
  int bond_dim_init{};
  int bond_dim_increase{};
  int bond_dim_max{};
  int sweep_bound{};
  std::vector<double> v_value{};
  std::vector<double> v_weight{};
  double auxi_height{};
  double reltol{};
  bool fullPiv{};
  int error_type{};
  int error_eval{};
  double convergence_bound{};
  int convergence_iter{};
  double decay_rate{};
};

struct simulation_params_t {
  debug_t debug{};
  size_t inch_start_index{};
  size_t inch_end_index{};
  double tau_max{};
  double tau_split{};
  double tau_split_ratio{};
  long n_tau_linear{};
  int order_Chebyshev{};
  long n_tot{};
  interpolation_type interp_type{};
  std::vector<double> grid{};
  std::vector<double> grid_linear{};
  std::vector<int> order_list {};
  std::vector<int> order_list_first {};
  int bl_index{};
  int subspace_index{};
  bool use_bare_propagator{};
  int eval_type{0}; // 0 for propagator/partition function, 1 for greens function
  int n_skip{0};
  std::vector<int> gf_index{};
};

struct VectorHash {
    std::size_t operator()(const std::vector<double>& v) const {
        std::size_t hash = v.size();
        for (const auto& i : v) {
            hash ^= std::hash<double>()(i) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};


struct eval_params_t {
  cv_func change_variable;
  jb_func jacobian; 
  std::unordered_map<int, std::vector<std::pair<std::vector<int>, std::vector<int>>>> phi_pair_order_cache;
  std::unordered_map<std::vector<double>, std::vector<std::pair<std::vector<int>, std::vector<int>>>,VectorHash> phi_pair_iota_cache;
};

struct statistics_per_inch_t {
  std::vector<long> func_evals_order {};
  std::vector<long> warning_same_time_order {};
  std::vector<long> warning_tau_split_order {};
  std::vector<long> warning_tau_max_order {};
  std::vector<double> max_diff_order {}; 
  std::vector<double> max_auxi_height_order {};
  std::vector<double> max_error_order {};
  std::vector<double> u_tau_sum_order {}; // this is the abs sum of all matrix element and all tau within the same inchworm/bare step
  std::vector<double> time_order {};
  std::vector<long> nTCI_order {};
  std::vector<double> integral_max_order {};
};

struct simulation_results_t {
  double Z_bath            = 0;
  double Z_bath_correction = 0;
  double Z_imp_correction  = 0;
  double partition_function{};
  u_tau_t u_tau{};
  g_tau_t G_tau{};
  interpolator_t<scalar_t> u_interpolator{};
  frame_t u_tau_max_zeroth_order{};
  double partition_function_ref{};
  u_tau_t u_tau_ref{};
  g_tau_t G_tau_ref{};
  interpolator_t<scalar_t> u_interpolator_ref{};
  frame_t u_tau_zeroth_order{};

  std::vector<std::vector<double>> integral_list = {};
  std::vector<statistics_per_inch_t> statistics{};
};

template <typename T> struct Loop {
  std::string name;
  T container;
  std::vector<double> value = {};
  Loop(std::string name, T container) : name(name), container(container) {}
  Loop &operator=(const Loop &other) {
    if (this != &other) { // protect against self-assignment
      this->name      = other.name;
      this->container = other.container;
    }
    return *this;
  }
};

class ModeBase {
  public:
  ModeBase() {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  }
  virtual void init(std::string json_file_path, std::string hyb_file_path) {
    NVTX_RANGE("init", 0);
    try {
      read_json_parameters(json_file_path);
    } catch (std::exception &e) {
      std::cerr << "Error in reading json file: " << e.what() << std::endl;
      std::exit(EXIT_FAILURE);
    }
    try {
      prepare_input(hyb_file_path);
    } catch (std::exception &e) {
      std::cerr << "Error in preparing input: " << e.what() << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  virtual void run() = 0;
  virtual void validate_input();
  virtual void print_summary();
  virtual void prepare_eval();
  virtual void evaluate(std::vector<std::vector<double>> const &unsummed_input         = std::vector<std::vector<double>>(),
                        bool is_first_interval = false, size_t inchworm_index = 0);
  virtual void evaluate_propagator()      = 0;
  virtual void evaluate_greens_function() = 0;
  virtual ~ModeBase() {}
  std::string mode_name{};
  int rank;
  int size;

  // friend function (saving files->save.hpp)
  friend void h5_save_params(const ModeBase *mode, h5::group h5group, std::string subgroup_name);
  friend void h5_save_propagator(const ModeBase *mode, h5::group h5group, std::string subgroup_name);
  friend void h5_save_cheb_coeff(const ModeBase *mode, h5::group h5group, std::string subgroup_name);
  friend void h5_save_propagator_ref(const ModeBase *mode, h5::group h5group, std::string subgroup_name);
  friend void h5_save_cheb_coeff_ref(const ModeBase *mode, h5::group h5group, std::string subgroup_name);
  friend void h5_save_gf(const ModeBase *mode, h5::group h5group, std::string subgroup_name, g_tau_t const &G_tau);
  friend void h5_save_statistics(const ModeBase *mode, h5::group h5group, std::string subgroup_name);

  protected:
  // parameters for all modes
  global_params_t gp{};
  constr_params_t cp{};
  model_params_t mp{};
  tci_params_t tp{};
  simulation_params_t sp{};
  simulation_results_t sr{};
  eval_params_t ep{};
  void read_json_parameters(std::string json_file_path);
  hyb_tau_t read_hyb_function(std::string hyb_file_path, model_params_t const &mp, constr_params_t const &cp);
  void prepare_input(std::string hyb_file_path);
  void clear_tci_results();
};

class ModeDebug : public ModeBase {
  public:
  ModeDebug() : ModeBase() { mode_name = "debug"; }
  void run() override;
  void validate_input() override;
  void evaluate_propagator() override;
  void evaluate_greens_function() override;
  void print_summary() override;
};

class ModeInchworm : public ModeBase {
  public:
  ModeInchworm() : ModeBase() { mode_name = "inchworm"; }
  void run() override;
  void validate_input() override;
  void evaluate_propagator() override;
  void evaluate_greens_function() override;
  void print_summary() override;
};

class ModeBare : public ModeBase {
  public:
  ModeBare() : ModeBase() { mode_name = "bare"; }
  void run() override;
  void validate_input() override;
  void evaluate_propagator() override;
  void evaluate_greens_function() override;
  void print_summary() override;
};
