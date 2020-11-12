#pragma once
#include "../qmc_config_data.hpp"
#include "../../container_set.hpp"

namespace inchworm::measures {

  /// g_frame measurement
  struct g_frame {

    // Constructor
    g_frame(params_t const &params, qmc_config_data_t const &qmc_config_data_, single_step_results_t &results_);

    // Invoke a single measurement
    void accumulate(scalar_t sign);

    // Finalize
    void collect_results(mpi::communicator const &comm);

    private:
    // The Monte-Carlo configuration
    qmc_config_data_t const &data;
    single_step_results_t &results;
  };

} // namespace inchworm::measures
