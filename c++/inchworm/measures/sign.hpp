#pragma once
#include "../qmc_data.hpp"
#include "../container_set.hpp"

namespace inchworm::measures {

  /// sign measurement
  struct sign {

    // Constructor
    // sign(params_t const &params, qmc_data_t &qmc_data_);
    sign(params_t const &params, qmc_data_t const &qmc_data_, container_set &results);

    // Invoke a single measurement
    void accumulate(mc_weight_t sign);

    // Finalize
    void collect_results(mpi::communicator const &comm);

    private:
    // The Monte-Carlo configuration
    qmc_data_t const &qmc_data;

    container_set &results;

    mc_weight_t average_sign_;

    long count;
  };

} // namespace inchworm::measures
