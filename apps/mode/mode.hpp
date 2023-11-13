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

class BaseMode {
  public:
  BaseMode() {}
  virtual void init(std::string jsonFilePath) {
    readJsonParameters(jsonFilePath, debug, cp, n_site, epsilon, theta, n_bath, n_spin, U, mu, t, tau_max, tau_split, n_GK, bond_dim, sweep_bound,
                       order_list, tci_prrlu, error_bound, bl_index, subspace_index);
  }
  void constructHubbard() {
    // prepare input
    std::tie(vi, wi) = selectQuadratureGK(n_GK, 0, 1);
    std::tie(Delta_tau, ad_imp, u_tau, G_tau) = test_setup(n_site, n_bath, n_spin, U, mu, t, cp, theta, epsilon);
    u_interpolator                            = interpolator_t<scalar_t>(u_tau, u_tau[0].mesh().size());
    u_tau_max_zeroth_order                    = u_interpolator(tau_max - tau_split) * u_interpolator(tau_split); //oder 0 result
    n_bl                                      = cp.gf_struct.size();
    all_d_ops.resize(n_bl, {});
    all_d_dag_ops.resize(n_bl, {});
    for (auto [bl, bl_pair] : enumerate(cp.gf_struct)) {
      auto [bl_name, bl_size] = bl_pair;
      block_shape.push_back(bl_size);
    }
    fops = fundamental_operator_set{cp.gf_struct};
    for (auto [bl, bl_pair] : enumerate(cp.gf_struct)) {
      auto [bl_name, bl_size] = bl_pair;
      all_d_ops[bl].clear();
      all_d_dag_ops[bl].clear();
      for (auto idx : range(bl_size)) {
        all_d_ops[bl].emplace_back(0.0, false, fops[{bl_name, idx}], bl, idx);
        all_d_dag_ops[bl].emplace_back(0.0, true, fops[{bl_name, idx}], bl, idx);
      }
    }
    if (debug) {
      std::cout << "Delta_tau shape:" << std::endl;
      printBlockShape(Delta_tau);
      std::cout << "G_tau shape:" << std::endl;
      printBlockShape(G_tau);
      std::cout << "u_tau shape:" << std::endl;
      printBlockShape(u_tau);
    }
  }

  virtual void runSingleElement() = 0;
  virtual ~BaseMode() {}

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
};

class ModeExplicitSum : public BaseMode {
  public:
  ModeExplicitSum() : BaseMode() {}
  void runSingleElement() override;
};

class ModeUseNormPivots : public BaseMode {
  public:
  ModeUseNormPivots() : BaseMode() {}
  void runSingleElement() override;
};

class ModeFullFactorization : public BaseMode {
  public:
  ModeFullFactorization() : BaseMode() {}
  void runSingleElement() override;
};

class ModeVertexFactorization : public BaseMode {
  public:
  ModeVertexFactorization() : BaseMode() {}
  void runSingleElement() override;
};

class ModeNestedTCI : public BaseMode {
  public:
  ModeNestedTCI() : BaseMode() {}
  void runSingleElement() override;
  void init(std::string jsonFilePath) override;

  private:
  // parameters for nested TCI only
  bool debug_iota{};
  bool tci_prrlu_iota{};
  int bond_dim_iota{};
  int sweep_bound_iota{};
  double error_bound_iota{};
};

class ModeReusePivots : public BaseMode {
  public:
  ModeReusePivots() : BaseMode() {}
  void runSingleElement() override;
};

class ModePartitionFactorization : public BaseMode {
  public:
  ModePartitionFactorization() : BaseMode() {}
  void runSingleElement() override;
};

class ModeCombineFactorization : public BaseMode {
  public:
  ModeCombineFactorization() : BaseMode() {}
  void runSingleElement() override;
};
