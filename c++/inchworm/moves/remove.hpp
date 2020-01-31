#pragma once
#include "./../qmc_config.hpp"

#include <triqs/mc_tools/random_generator.hpp>
#include <inchworm/diagram/inclusion_exclusion.hpp>

namespace inchworm::moves {

  /// A simple Monte-Carlo move
  struct remove {

    /// The Monte-Carlo configuration
    qmc_config_t &data;

    /// The random number generator
    triqs::mc_tools::random_generator &rng;

    /// Attempt vertex removing
    mc_weight_t attempt();

    /// Accept vertex removing
    mc_weight_t accept();

    /// Reject vertex removing
    void reject();

    /// w_hyb of proposed remove
    hybridization_scalar_t new_w_hyb = 1.0;

    /// w_loc of proposed remove
    double new_w_loc = 1.0;
  };

} // namespace inchworm::moves
