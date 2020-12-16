#pragma once

#include "../params.hpp"
#include "../qmc_data.hpp"
#include "../container_set.hpp" // qmc_results_t

namespace inchworm::measures {

  /// Measure of the average perturbation order
  struct average_order {

    average_order(params_t const &, qmc_data_t const &qmc_data_, qmc_results_t &results);

    /// Accumulate average sign
    void accumulate(scalar_t);

    /// Reduce and normalize
    void collect_results(mpi::communicator const &comm);

    private:
    // The Monte-Carlo configuration and data
    qmc_data_t const &qmc_data;

    // Reference to double for accumulation
    double &average_order_ref;

    // Accumulation counter
    long long N = 0;
  };

} // namespace inchworm::measures
