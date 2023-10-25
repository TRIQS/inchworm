#pragma once

#include "../params.hpp"
#include "../config.hpp"
#include "../container_set.hpp" // qmc_results_t

namespace inchworm::measures {

  // Measure the histogram of perturbation order
  struct order_histogram {

    order_histogram(params_t const &, config_t const &config, qmc_results_t &results);

    /// Accumulate perturbation order into histogram
    void accumulate(scalar_t);

    /// Reduce and normalize
    void collect_results(mpi::communicator const &comm);

    private:
    // The Monte-Carlo configuration and data
    config_t const &config;

    // Reference to accumulation vector
    std::vector<double> &order_histogram_ref;

    // Accumulation counter
    long long N = 0;
  };

} // namespace inchworm::measures
