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

    // Errors of the frame_[k][0](0,0) components
    std::vector<scalar_t>& err_frame_by_order;

    // Reference to the accumulation vector
    std::vector<frame_t> &acc_frame_by_order;

    // A zero initialized frame
    frame_t zero_frame;

    // The vector of accumulators for the error analysis
    std::vector<accumulator<scalar_t>> lin_acc_by_order;
    accumulator<scalar_t> empty_lin_acc = {0.0, 1000, -1};
  };

} // namespace inchworm::measures
