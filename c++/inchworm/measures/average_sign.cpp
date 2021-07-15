#include "./average_sign.hpp"

namespace inchworm::measures {

  average_sign::average_sign(params_t const &, config_t const &config, qmc_results_t &results)
     : config(config), average_sign_ref(results.average_sign) {
    average_sign_ref = 0.0;
  }

  void average_sign::accumulate(scalar_t sign) {
    average_sign_ref += sign;
    ++N;
  }

  void average_sign::collect_results(mpi::communicator const &comm) {
    N = mpi::all_reduce(N, comm);

    // Reduce and normalize
    average_sign_ref = mpi::all_reduce(average_sign_ref, comm);
    average_sign_ref = average_sign_ref / N;
  }

} // namespace inchworm::measures
