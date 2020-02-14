#include "./average_k.hpp"

namespace inchworm::measures {

  average_k::average_k(params_t const &, qmc_config_data_t const &qmc_config_data_, single_step_results_t &results_)
     : qmc_config_data(qmc_config_data_), results(results_) {}

  void average_k::accumulate(scalar_t sign) {
    ave_k += qmc_config_data.config.size();
    ++count;
  }

  void average_k::collect_results(mpi::communicator const &comm) {
    ave_k = mpi::all_reduce(ave_k, comm);
    count = mpi::all_reduce(count, comm);
    ave_k /= count;

    results.average_k = ave_k;
  }

} // namespace inchworm::measures
