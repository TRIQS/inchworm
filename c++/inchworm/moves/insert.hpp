#pragma once
#include "./../qmc_config.hpp"

#include <triqs/mc_tools/random_generator.hpp>
#include <inchworm/diagram/inclusion_exclusion.hpp>

namespace inchworm::moves {

  /// A simple Monte-Carlo move
  struct insert {

    /// The Monte-Carlo configuration
    qmc_config_t &data;

    /// The random number generator
    triqs::mc_tools::random_generator &rng;

    /// Attempt vertex insertion
    mc_weight_t attempt();

    /// Accept vertex insertion
    mc_weight_t accept();

    /// Reject vertex insertion
    void reject();

    /// w_hyb of proposed insert
    hybridization_scalar_t new_w_hyb = 1.0;

    /// w_loc of proposed insert
    double new_w_loc = 1.0;
  };

} // namespace inchworm::moves
