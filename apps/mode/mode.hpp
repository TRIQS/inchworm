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

class base_mode {
  public:
  base_mode() {}
  virtual void init(std::string json_file_path) {
    read_json_parameters(json_file_path, debug, cp, n_site, epsilon, theta, n_bath, n_spin, U, mu, t, tau_max, tau_split, n_GK, bond_dim, sweep_bound,
                       order_list, tci_prrlu, error_bound, bl_index, subspace_index);
  }
  void construct_Hubbard();
  void print_summary();
  virtual void run_single_element() = 0;
  virtual ~base_mode() {}

  protected:
  // input parameters
  bool debug{};
  constr_params_t cp{};
  mat_t theta;
  vec_t epsilon;
  int n_site{};
  int n_bath{};
  int n_spin{};
  int n_GK{};
  int bond_dim{};
  int sweep_bound{};
  double U{};
  double mu{};
  double t{};
  double tau_max{};
  double tau_split{};
  std::vector<int> order_list = {};
  bool tci_prrlu{};
  double error_bound{};
  int bl_index{};
  int subspace_index{};
  // constructed initial data
  std::vector<double> vi{};
  std::vector<double> wi{};
  hyb_tau_t Delta_tau{};
  atom_diag ad_imp{};
  u_tau_t u_tau{};
  g_tau_t G_tau{};
  interpolator_t<scalar_t> u_interpolator{};
  frame_t u_tau_max_zeroth_order{};
  std::vector<int> block_shape{};
  std::vector<std::vector<fop_t>> all_d_ops{};
  std::vector<std::vector<fop_t>> all_d_dag_ops{};
  long n_bl{};
  fundamental_operator_set fops{};
  //results
  std::vector<double> integral_order_list   = {};
  std::vector<double> calculation_time_list = {};
};

class ModeExplicitSum : public base_mode {
  public:
  ModeExplicitSum() : base_mode() {}
  void run_single_element() override;
};

class ModeUseNormPivots : public base_mode {
  public:
  ModeUseNormPivots() : base_mode() {}
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

class ModeNestedTCI : public base_mode {
  public:
  ModeNestedTCI() : base_mode() {}
  void run_single_element() override;
  void init(std::string json_file_path) override;

  private:
  // parameters for nested TCI only
  bool debug_iota{};
  bool tci_prrlu_iota{};
  int bond_dim_iota{};
  int sweep_bound_iota{};
  double error_bound_iota{};
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
