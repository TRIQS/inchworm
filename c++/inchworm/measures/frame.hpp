#pragma once

#include "../params.hpp"
#include "../config.hpp"
#include "../container_set.hpp" // qmc_results_t

namespace inchworm::measures {

  /// Measurement of a single frame of a propagator or Green function
  struct frame {

    frame(params_t const &, config_t const &config, frame_t const &frame, qmc_results_t &results);

    /// Invoke a single measurement
    void accumulate(scalar_t);

    /// Finalize
    void collect_results(mpi::communicator const &comm);

    private:
    // Print verbosity
    int verbosity;

    // The Monte-Carlo configuration
    config_t const &config;

    // The current frame
    frame_t const &frame_;

    // Errors of the (0,0) component for each frame[bl]
    std::vector<scalar_t> &errs_frame;

    // References to the accumulation frames
    frame_t &acc_frame;
    frame_t &acc_frame_zeroth_order;

    // The scalar accumulator for the auto-correlation analysis
    accumulator<scalar_t> log_acc = {0.0, -1, 0};

    // The scalar accumulator for the error analysis
    std::vector<accumulator<scalar_t>> lin_acc;

    // The number of samples
    long long N_samples = 0;
  };

} // namespace inchworm::measures
