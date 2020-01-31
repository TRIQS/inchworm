#pragma once
#include "../qmc_config.hpp"
#include "../container_set.hpp"

namespace inchworm::measures {

  /// U_frame measurement
  struct U_frame {

    // Constructor
    sign(params_t const &params, qmc_config_t const &qmc_config_);

    // Invoke a single measurement
    void accumulate(mc_weight_t sign);

    // Finalize
    void collect_results(mpi::communicator const &comm);

    private:
    // The Monte-Carlo configuration
    qmc_config_t const &qmc_config;

    mc_weight_t average_sign_; 
    propagator_frame U_frame_;
  };

} // namespace inchworm::measures
