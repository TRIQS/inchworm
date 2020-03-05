#pragma once
#include "../qmc_config_data.hpp"
#include "../../diagram/print.hpp"
#include "../../u_frame.hpp"

#include <triqs/mc_tools/random_generator.hpp>
#include <inchworm/diagram/inclusion_exclusion.hpp>

namespace inchworm::moves {

  /// A simple Monte-Carlo move
  struct insert {

    /// Attempt vertex insertion
    scalar_t attempt();

    /// Accept vertex insertion
    scalar_t accept();

    /// Reject vertex insertion
    void reject() {}

    /// Constructor
    insert(qmc_config_data_t &data, qmc_params_t const &qmc_params, triqs::mc_tools::random_generator &rng) : data(data), params(qmc_params), rng(rng) {}

    private:
    /// The Monte-Carlo configuration
    qmc_config_data_t &data;

    /// The Monte-Carlo parameters
    qmc_params_t const &params;

    /// The random number generator
    triqs::mc_tools::random_generator &rng;

    /// c/cdag lists of proposed insert
    config_t proposed_config; // proposed configuration of c and cdag

    /// weights of the proposed configuration
    weights_t proposed_w;

    /// container of the calculated time frame of the propagator
    u_frame_t proposed_u_frame = make_zero_propagator_frame(params.h_diag);

    int proposed_sign;
  };

} // namespace inchworm::moves
