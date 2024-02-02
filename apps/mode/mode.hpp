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
};

struct tci_params_t {
  int n_GK{};
  int mapping_v{}; //0 for the mapping in Phys. Rev. B 107, 245135, 1 for the mapping in 	arXiv:2310.16957
  bool tci_prrlu{};
  int bond_dim{};
  int sweep_bound{};
  std::vector<double> v_value{};
  std::vector<double> v_weight{};
  double auxi_height{};
  double reltol{};
  double integral_lower_bound{};
  double convergence_bound{};
  int convergence_iter{};
};

struct simulation_params_t {
  debug_t debug{};
  double tau_max{};
  double tau_split{};
  std::vector<int> order_list = {};
  int bl_index{};
  int subspace_index{};
};

struct simulation_results_t {
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

  std::vector<double> integral_list         = {};
  std::vector<double> calculation_time_list = {};
  std::vector<double> pretrain_time_list    = {};
  std::vector<double> find_pivot_time_list  = {};
  std::vector<double> train_time_list       = {};
};

class ModeBase {
  public:
  ModeBase() {}
  virtual void init(std::string json_file_path) {
    read_json_parameters(json_file_path);
    prepare_input();
  }
  virtual void print_summary();
  virtual void run() = 0;
  virtual ~ModeBase() {}
  std::string mode_name{};

  protected:
  // parameters for all modes
  global_params_t gp{};
  constr_params_t cp{};
  model_params_t mp{};
  tci_params_t tp{};
  simulation_params_t sp{};
  simulation_results_t sr{};
  void read_json_parameters(std::string json_file_path);
  void prepare_input();
};

class ModeDebug : public ModeBase {
  public:
  ModeDebug() : ModeBase() { mode_name = "debug"; }
  void run() override;
};

class ModeInchworm : public ModeBase {
  public:
  ModeInchworm() : ModeBase() { mode_name = "inchworm"; }
  void run() override;
};

class ModeBare : public ModeBase {
  public:
  ModeBare() : ModeBase() { mode_name = "bare"; }
  void run() override;
};
