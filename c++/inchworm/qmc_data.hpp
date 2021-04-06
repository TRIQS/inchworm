#pragma once

#include "types.hpp"
#include "u_frame.hpp"

namespace inchworm {

  struct config_t {
    std::vector<fop_t> d_list, d_dag_list;                    // list of d/d_dag time ordered
    std::vector<std::vector<fop_t>> d_bl_list, d_dag_bl_list; // list of d/d_dag by block, insertion ordered

    config_t(long n_bl) : d_bl_list(n_bl), d_dag_bl_list(n_bl) {}

    long size() const { return d_list.size(); }
    long size(long bl) const { return d_bl_list[bl].size(); }

    bool try_insert(fop_t const &ddag, fop_t const &d);
    bool try_erase(long bl, long i_dag, long i);
    bool try_double_insert(fop_t const &d_dag1, fop_t const &d1, fop_t const &d_dag2, fop_t const &d2);
    bool try_double_erase(long bl1, long i1_dag, long i1, long bl2, long i2_dag, long i2);
  };

  struct weights_t {
    scalar_t imp; // Impurity weight - Frobenius norm of the current (propagator or green function) frame
    scalar_t hyb; // Hybridization weight - Value of the determinant in cthyb or its equivalent for the inchworm
  };

  /// The Monte-Carlo Configuration structure
  struct qmc_data_t {
    // last accepted paraemeters
    config_t config;   // last accepted configuration of d and d_dag
    weights_t weights; // weight values of the last accepted configuration
    int sign = 1;      // sign of the last accepted configuration
    frame_t frame;     // the configuration (propagator or green function) frame

    qmc_data_t(gf_struct_t const &gf_struct, frame_t const &frame) : config(gf_struct.size()), weights{frobenius_norm(frame), 1.0}, frame{frame} {}
  };

  // static parameters of the Monte Carlo simulation
  struct qmc_params_t {
    double tau_max;               // similar to beta, but configuration here does not always goes up to beta. 0 < tau_max <= beta
    double tau_split;             // in the inchworm, this should be the tau_max of the previous inching. 0 < tau_split <= tau_max
    bool use_bare_propagator;     // True only for the first iteration of the inchworm calculation
    MODE mode;                    // The sampling mode, either PROPAGATOR or GREENFUNCTION
    std::optional<int> max_order; // The maximum perturbation order [optional]
  };
} // namespace inchworm
