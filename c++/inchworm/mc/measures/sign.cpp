#include "./sign.hpp"

namespace inchworm::measures {

  sign::sign(params_t const &, qmc_config_data_t const &qmc_config_data_, single_step_results_t &results_)
     : qmc_config_data(qmc_config_data_), results(results_) {}

  void sign::accumulate(scalar_t sign) {
    average_sign += sign;
    ++count;
  }

  void sign::collect_results(mpi::communicator const &comm) {
    average_sign = mpi::all_reduce(average_sign, comm);
    count        = mpi::all_reduce(count, comm);
    average_sign /= count;

    results.average_sign = average_sign;
  }

} // namespace inchworm::measures
