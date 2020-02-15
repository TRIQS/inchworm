#pragma once
#include "../qmc_config_data.hpp"
#include "../../container_set.hpp"

namespace inchworm::measures {

  /// sign measurement
  struct sign {

    // Constructor
    // sign(params_t const &params, qmc_config_data_t &qmc_config_data_);
    sign(params_t const &, qmc_config_data_t const &qmc_config_data_, single_step_results_t &results_);

    // Invoke a single measurement
    void accumulate(scalar_t sign);

    // Finalize
    void collect_results(mpi::communicator const &comm);

    private:
    // The Monte-Carlo configuration
    qmc_config_data_t const &qmc_config_data;

    single_step_results_t &results;

    scalar_t average_sign = 0;

    scalar_t z = 0;
  };

} // namespace inchworm::measures
