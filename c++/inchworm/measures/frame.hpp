#pragma once

#include "../params.hpp"
#include "../qmc_data.hpp"
#include "../container_set.hpp" // qmc_results_t

namespace inchworm::measures {

  /// Measurement of a single frame of a propagator or Green function
  struct frame {

    frame(params_t const &, qmc_data_t const &qmc_data_, qmc_results_t &results);

    /// Invoke a single measurement
    void accumulate(scalar_t);

    /// Finalize
    void collect_results(mpi::communicator const &comm);

    private:
    // Print verbosity
    int verbosity;

    // The Monte-Carlo configuration and data
    qmc_data_t const &qmc_data;

    // References to the accumulation frames
    frame_t &frame_ref;
    frame_t &frame_0th_order_ref;

    // The scalar accumulator for the auto-correlation analysis
    accumulator<scalar_t> acc = {0.0, -1};
  };

} // namespace inchworm::measures
