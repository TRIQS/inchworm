#pragma once
#include "../qmc_config.hpp"
#include "../container_set.hpp"

namespace inchworm::measures {

  /// A simple measurement
  struct simple {

    // Constructor
    simple(params_t const &params, qmc_config_t const &qmc_config_, container_set &results);

    // Invoke a single measurement
    void accumulate(mc_weight_t sign);

    // Finalize
    void collect_results(mpi::communicator const &comm);

    private:
    // The Monte-Carlo configuration
    qmc_config_t const &qmc_config;
  };

} // namespace inchworm::measures
