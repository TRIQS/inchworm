#pragma once

#include "../params.hpp"
#include "../qmc_data.hpp"
#include "../container_set.hpp" // qmc_results_t

namespace inchworm::measures {

  // Measure the histogram of perturbation order
  struct order_histogram {

    order_histogram(params_t const &, qmc_data_t const &qmc_data_, qmc_results_t &results);

    /// Accumulate perturbation order into histogram
    void accumulate(scalar_t);

    /// Reduce and normalize
    void collect_results(mpi::communicator const &comm);

    private:
    // The Monte-Carlo configuration and data
    qmc_data_t const &qmc_data;

    // Reference to accumulation vector
    std::vector<double> &order_histogram_ref;

    // Accumulation counter
    long long N = 0;
  };

} // namespace inchworm::measures
