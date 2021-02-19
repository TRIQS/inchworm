#pragma once

#include "../params.hpp"
#include "../qmc_data.hpp"
#include "../container_set.hpp" // qmc_results_t

namespace inchworm::measures {

  /// Measurement of a autocorr based on the trace of the current frame
  struct autocorr {

    autocorr(params_t const &, qmc_data_t const &qmc_data_, qmc_results_t &results);

    /// Invoke a single measurement
    void accumulate(scalar_t);

    /// Finalize
    void collect_results(mpi::communicator const &comm);

    private:
    // Print verbosity
    int verbosity;

    // The Monte-Carlo configuration and data
    qmc_data_t const &qmc_data;

    // The Result Container
    qmc_results_t &results;

    // The scalar log-bin accumulator for the auto-correlation analysis
    accumulator<scalar_t> acc = {0.0, -1};
  };

} // namespace inchworm::measures
