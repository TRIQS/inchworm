#pragma once
#include "../qmc_config.hpp"
#include "../container_set.hpp"

namespace inchworm::measures {

  /// U_frame measurement
  struct U_frame {

    // Constructor
    U_frame(params_t const &params, qmc_config_t const &qmc_config_, container_set &results_);

    // Invoke a single measurement
    void accumulate(mc_weight_t sign);

    // Finalize
    void collect_results(mpi::communicator const &comm);

    private:
    // The Monte-Carlo configuration
    qmc_config_t const &qmc_config;
    container_set &results;

    mc_weight_t average_sign_; 
    propagator_frame U_frame_;
  };

} // namespace inchworm::measures
