#pragma once
#include "../qmc_config_data.hpp"
#include "../../diagram/print.hpp"
#include "../../u_frame.hpp"

#include <triqs/mc_tools/random_generator.hpp>
#include <inchworm/diagram/inclusion_exclusion.hpp>
#include <inchworm/diagram/proper_enum.hpp>

namespace inchworm::moves {

  /// A simple Monte-Carlo move
  struct remove {

    /// Attempt vertex removing
    scalar_t attempt();

    /// Accept vertex removing
    scalar_t accept();

    /// Reject vertex removing
    void reject() {}

    /// Constructor:
    remove(qmc_config_data_t &data, params_t const &params, qmc_params_t const &qmc_params, triqs::mc_tools::random_generator &rng)
       : data(data), params(qmc_params), rng(rng) {
      proposed_g_frame = make_zero_green_frame(params.gf_struct);
    }

    private:
    /// The Monte-Carlo configuration
    qmc_config_data_t &data;

    /// The Monte-Carlo parameters
    qmc_params_t const &params;

    /// The random number generator
    triqs::mc_tools::random_generator &rng;

    /// c/d_dag lists of proposed remove
    config_t proposed_config; // proposed configuration of c and d_dag

    /// weights of the proposed configuration
    weights_t proposed_w;

    /// container of the calculated time frame of the propagator
    u_partial_t proposed_u_partial;

    /// green function container to accumulate into
    g_frame_t proposed_g_frame;

    int proposed_sign;
  };

} // namespace inchworm::moves
