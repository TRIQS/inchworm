#pragma once
#include "../qmc_data.hpp"
#include "../container_set.hpp"

namespace inchworm::measures {

  /// u_frame measurement
  struct u_frame {

    // Constructor
    u_frame(params_t const &params, qmc_data_t const &qmc_data_, container_set &results_);

    // Invoke a single measurement
    void accumulate(scalar_t sign);

    // Finalize
    void collect_results(mpi::communicator const &comm);

    private:
    // The Monte-Carlo configuration
    qmc_data_t const &qmc_data;
    container_set &results;

    scalar_t average_sign_; 
    u_frame_t u_frame_;
  };

} // namespace inchworm::measures
