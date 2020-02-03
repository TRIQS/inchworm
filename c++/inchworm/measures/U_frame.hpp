#pragma once
#include "../qmc_data.hpp"
#include "../container_set.hpp"

namespace inchworm::measures {

  /// U_frame measurement
  struct U_frame {

    // Constructor
    U_frame(params_t const &params, qmc_data_t const &qmc_data_, container_set &results_);

    // Invoke a single measurement
    void accumulate(mc_weight_t sign);

    // Finalize
    void collect_results(mpi::communicator const &comm);

    private:
    // The Monte-Carlo configuration
    qmc_data_t const &qmc_data;
    container_set &results;

    mc_weight_t average_sign_; 
    u_frame_t U_frame_;
  };

} // namespace inchworm::measures
