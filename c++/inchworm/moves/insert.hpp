#pragma once
#include "./../qmc_config.hpp"

#include <triqs/mc_tools/random_generator.hpp>
#include <inchworm/diagram/inclusion_exclusion.hpp>

namespace inchworm::moves {

  /// A simple Monte-Carlo move
  struct insert {

    /// w_hyb of proposed insert
    propagator_frame new_U_frame;

    /// w_hyb of proposed insert
    hybridization_scalar_t new_w_hyb;
    
    /// w_loc of proposed insert
    double new_w_loc;

    /// The Monte-Carlo configuration
    qmc_config_t &qmc_config;

    /// The random number generator
    triqs::mc_tools::random_generator &rng;

    /// Attempt vertex insertion
    mc_weight_t attempt();

    /// Accept vertex insertion
    mc_weight_t accept();

    /// Reject vertex insertion
    void reject();
  };

} // namespace inchworm::moves
