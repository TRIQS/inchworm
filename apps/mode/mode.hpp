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
#include "../hubbard.hpp"
#include "../utility.hpp"

struct global_params_t {
  std::string target{};
  std::string integrand{};
  std::string integral_variable{};
  std::string tci_shape{};
  std::string trick{};
  int model_type{}; // 0 for discrete bath, 1 for continuous bath, 2 for bethe lattice
  bool do_segment = false;
  std::string output_prefix{};
  bool exact_sum                   = false;
  int unsummed_tci                = 0; // 0: all indices are summed; 1: the first index is not summed; 2: the first two indices are not summed ... 
  double energy_shift              = 0.0;
  double Z_energy_shift_correction = 1.0;
};

struct model_params_t {
  int n_site{};
  int n_bath{};
  int n_spin{};
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
  int n_omega_bethe{};
};

struct tci_params_t {
  int n_GK{};
  int mapping_v{}; //0 for the mapping in Phys. Rev. B 107, 245135, 1 for the mapping in 	arXiv:2310.16957
  int tci_prrlu{};
  int bond_dim{};
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
  double integral_lower_bound{};
  double decay_rate{};
};

struct simulation_params_t {
  debug_t debug{};
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
  frame_t u_tau_zeroth_order_ref{};
  frame_t u_tau_zeroth_order_bare{};
  frame_t u_tau_zeroth_order{};
  double partition_function_zeroth_order_ref{};

  std::vector<std::vector<double>> integral_list = {};
  std::vector<double> calculation_time_list      = {};
  std::vector<double> pretrain_time_list         = {};
  std::vector<double> find_pivot_time_list       = {};
  std::vector<double> train_time_list            = {};
  double total_time=0;
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
  ModeBase() {}
  virtual void init(std::string json_file_path, std::string hyb_file_path) {
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
  virtual void print_summary();
  virtual void run() = 0;
  virtual void validate_input();
  virtual void evaluate(std::vector<std::vector<double>> const &unsummed_input         = std::vector<std::vector<double>>(),
                        std::vector<std::vector<std::vector<double>>> const &all_input = std::vector<std::vector<std::vector<double>>>(),
                        std::vector<std::vector<double>> const &all_weight             = std::vector<std::vector<double>>(),
                        bool is_first_interval = false);
  virtual void evaluate_propagator()      = 0;
  virtual void evaluate_greens_function() = 0;
  virtual ~ModeBase() {}
  std::string mode_name{};

  // friend function (saving files->save.hpp)
  friend void h5_save_params(const ModeBase *mode, h5::group h5group, std::string subgroup_name);
  friend void h5_save_propagator(const ModeBase *mode, h5::group h5group, std::string subgroup_name);
  friend void h5_save_cheb_coeff(const ModeBase *mode, h5::group h5group, std::string subgroup_name);
  friend void h5_save_propagator_ref(const ModeBase *mode, h5::group h5group, std::string subgroup_name);
  friend void h5_save_cheb_coeff_ref(const ModeBase *mode, h5::group h5group, std::string subgroup_name);
  friend void h5_save_gf(const ModeBase *mode, h5::group h5group, std::string subgroup_name, g_tau_t const &G_tau);

  protected:
  // parameters for all modes
  global_params_t gp{};
  constr_params_t cp{};
  model_params_t mp{};
  tci_params_t tp{};
  simulation_params_t sp{};
  simulation_results_t sr{};
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
};

class ModeInchworm : public ModeBase {
  public:
  ModeInchworm() : ModeBase() { mode_name = "inchworm"; }
  void run() override;
  void validate_input() override;
  void evaluate_propagator() override;
  void evaluate_greens_function() override;
};

class ModeBare : public ModeBase {
  public:
  ModeBare() : ModeBase() { mode_name = "bare"; }
  void run() override;
  void validate_input() override;
  void evaluate_propagator() override;
  void evaluate_greens_function() override;
};
