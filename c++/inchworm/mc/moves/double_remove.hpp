#pragma once
#include "../qmc_config_data.hpp"
#include "../../diagram/print.hpp"
#include "../../u_frame.hpp"

#include <triqs/mc_tools/random_generator.hpp>
#include <inchworm/diagram/inclusion_exclusion.hpp>
#include <inchworm/diagram/proper_enum.hpp>

namespace inchworm::moves {

  /// A simple Monte-Carlo move
  struct double_remove {

    /// Attempt vertex double removal
    scalar_t attempt();

    /// Accept vertex double removal
    scalar_t accept();

    /// Reject vertex double removal
    void reject() {}

    /// Constructor:
    double_remove(qmc_config_data_t &data, params_t const &params, qmc_params_t const &qmc_params, triqs::mc_tools::random_generator &rng)
       : data(data), params(qmc_params), rng(rng), gf_struct(params.gf_struct) {
      prop_g_frame = make_frame(params.gf_struct);
    }

    private:
    /// The Monte-Carlo configuration
    qmc_config_data_t &data;

    /// The Monte-Carlo parameters
    qmc_params_t const &params;

    /// The random number generator
    triqs::mc_tools::random_generator &rng;

    /// d/d_dag lists of proposed double_remove
    config_t prop_config; // proposed configuration of d and d_dag

    /// weights of the proposed configuration
    weights_t prop_weights;

    /// container of the calculated time frame of the propagator
    u_partial_t prop_u_partial;

    /// the green function structure
    gf_struct_t const &gf_struct;

    /// green function container to accumulate into
    frame_t prop_g_frame;

    int prop_sign;
  };

} // namespace inchworm::moves
