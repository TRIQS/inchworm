#pragma once

#include "../params.hpp"
#include "../config.hpp"
#include "../container_set.hpp" // qmc_results_t

namespace inchworm::measures {

  /// Measurement of a autocorr based on the trace of the current frame
  struct autocorr {

    autocorr(params_t const &, config_t const &config, frame_t const &frame, qmc_results_t &results);

    /// Invoke a single measurement
    void accumulate(scalar_t);

    /// Finalize
    void collect_results(mpi::communicator const &comm);

    private:
    // Print verbosity
    int verbosity;

    // The Monte-Carlo configuration and data
    config_t const &config;

    // The current frame
    frame_t const &curr_frame;

    // The Result Container
    qmc_results_t &results;

    // The scalar log-bin accumulator for the auto-correlation analysis
    accumulator<scalar_t> log_acc = {0.0, -1, 0};
  };

} // namespace inchworm::measures
