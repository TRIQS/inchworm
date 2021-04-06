#pragma once

#include "../params.hpp"
#include "../config.hpp"
#include "../container_set.hpp" // qmc_results_t

namespace inchworm::measures {

  /// For each diagram order, measure the frame of a propagator or Green function
  struct frame_by_order {

    frame_by_order(params_t const &params, config_t const &config, frame_t const &frame, qmc_results_t &results);

    // Invoke a single measurement
    void accumulate(scalar_t sign);

    // Finalize
    void collect_results(mpi::communicator const &comm);

    private:
    // The Monte-Carlo configuration and data
    config_t const &config;

    // The current frame
    frame_t const &frame_;

    // Reference to the accumulation vector
    std::vector<frame_t> &acc_frame_by_order;

    // A zero initialized frame
    frame_t zero_frame;
  };

} // namespace inchworm::measures
