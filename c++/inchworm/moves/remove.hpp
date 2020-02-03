#pragma once
#include "../qmc_data.hpp"

#include <triqs/mc_tools/random_generator.hpp>
#include <inchworm/diagram/inclusion_exclusion.hpp>

namespace inchworm::moves {

  /// A simple Monte-Carlo move
  struct remove {

    /// The Monte-Carlo configuration
    qmc_data_t &data;

    /// The random number generator
    triqs::mc_tools::random_generator &rng;

    /// Attempt vertex removing
    mc_weight_t attempt();

    /// Accept vertex removing
    mc_weight_t accept();

    /// Reject vertex removing
    void reject() {}

    /// Constructor:
    remove(qmc_data_t &data, triqs::mc_tools::random_generator &rng) : data(data), rng(rng), new_u_frame(data.h_diag) {}

    private:
    /// w_hyb of proposed remove
    hyb_scalar_t new_w_hyb = 1.0;

    /// w_loc of proposed remove
    double new_w_loc = 1.0;

    /// c/cdag lists of proposed remove
    config_t new_config; // proposed configuration of c and cdag

    /// weights of the proposed configuration
    weights_t new_w;

    /// container of the calculated time frame of the propagator
    u_frame_t new_u_frame;

  };

} // namespace inchworm::moves
