#pragma once
#include "../qmc_config.hpp"

namespace inchworm::measures {

  /// sign measurement
  struct sign {

    // Constructor
    sign(params_t const &params, qmc_config_t &qmc_config_);

    // Invoke a single measurement
    void accumulate(mc_weight_t sign);

    // Finalize
    void collect_results(mpi::communicator const &comm);

    private:
    // The Monte-Carlo configuration
    qmc_config_t const &qmc_config;

    mc_weight_t average_sign_; 

    long count;
  };

} // namespace inchworm::measures
