#pragma once

#include "../params.hpp"
#include "../config.hpp"
#include "../container_set.hpp" // qmc_results_t

namespace inchworm::measures {

  /// Measure of the average sign
  struct average_sign {

    average_sign(params_t const &, config_t const &, qmc_results_t &results);

    /// Accumulate average sign
    void accumulate(scalar_t);

    /// Reduce and normalize
    void collect_results(mpi::communicator const &comm);

    private:
    // Reference to double for accumulation
    scalar_t &average_sign_ref;

    // Accumulation counter
    long long N = 0;
  };

} // namespace inchworm::measures
