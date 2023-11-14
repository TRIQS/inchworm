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

enum debug_t {
  none, //0, no debug
  low,  //1, only important information is printed
  high  //2, all information is printed
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
};

struct tci_params_t {
  int n_GK{};
  bool tci_prrlu{};
  int bond_dim{};
  int sweep_bound{};
  double error_bound{};
  std::vector<double> vi{};
  std::vector<double> wi_v{};
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
  u_tau_t u_tau{};
  g_tau_t G_tau{};
  interpolator_t<scalar_t> u_interpolator{};
  frame_t u_tau_max_zeroth_order{};
  std::vector<double> integral_order_list   = {};
  std::vector<double> calculation_time_list = {};
};

class base_mode {
  public:
  base_mode() {}
  virtual void init(std::string json_file_path) { base_mode::read_json_parameters(json_file_path); }
  void prepare_input();
  void print_summary();
  virtual void run_single_element() = 0;
  virtual void read_json_parameters(std::string json_file_path);
  virtual ~base_mode() {}

  protected:
  // parameters for all modes
  constr_params_t cp{};
  model_params_t mp{};
  tci_params_t tp{};
  simulation_params_t sp{};
  simulation_results_t sr{};
};

class ModeExplicitSum : public base_mode {
  public:
  ModeExplicitSum() : base_mode() {}
  void run_single_element() override;
};

class ModeFullFactorization : public base_mode {
  public:
  ModeFullFactorization() : base_mode() {}
  void run_single_element() override;
};

class ModeVertexFactorization : public base_mode {
  public:
  ModeVertexFactorization() : base_mode() {}
  void run_single_element() override;
};

class ModeReusePivots : public base_mode {
  public:
  ModeReusePivots() : base_mode() {}
  void run_single_element() override;
};

class ModePartitionFactorization : public base_mode {
  public:
  ModePartitionFactorization() : base_mode() {}
  void run_single_element() override;
};

class ModeCombineFactorization : public base_mode {
  public:
  ModeCombineFactorization() : base_mode() {}
  void run_single_element() override;
};

class ModeNestedTCI : public base_mode {
  public:
  ModeNestedTCI() : base_mode() {}
  void run_single_element() override;
  void init(std::string json_file_path) override;
  void read_json_parameters(std::string json_file_path) override;

  private:
  // parameters for nested TCI only
  tci_params_t tp_iota{};
};
