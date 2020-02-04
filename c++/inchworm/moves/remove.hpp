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
    scalar_t attempt();

    /// Accept vertex removing
    scalar_t accept();

    /// Reject vertex removing
    void reject() {}

    /// Constructor:
    //remove(qmc_data_t &data, triqs::mc_tools::random_generator &rng) : data(data), rng(rng) { proposed_u_frame = init_propagator_frame(data.h_diag); }

    private:
    /// c/cdag lists of proposed remove
    config_t proposed_config; // proposed configuration of c and cdag

    /// weights of the proposed configuration
    weights_t proposed_w;

    /// container of the calculated time frame of the propagator
    u_frame_t proposed_u_frame = make_zero_propagator_frame(data.h_diag);

  };

} // namespace inchworm::moves
