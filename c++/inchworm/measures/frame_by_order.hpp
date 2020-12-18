#pragma once

#include "../params.hpp"
#include "../qmc_data.hpp"
#include "../container_set.hpp" // qmc_results_t

namespace inchworm::measures {

  /// For each diagram order, measure the frame of a propagator or Green function
  struct frame_by_order {

    frame_by_order(params_t const &params, qmc_data_t const &qmc_data_, qmc_results_t &results_);

    // Invoke a single measurement
    void accumulate(scalar_t sign);

    // Finalize
    void collect_results(mpi::communicator const &comm);

    private:
    // The Monte-Carlo configuration and data
    qmc_data_t const &qmc_data;

    // Reference to the accumulation vector
    std::vector<frame_t> &frame_by_order_ref;

    // A zero initialized frame
    frame_t zero_frame;
  };

} // namespace inchworm::measures
