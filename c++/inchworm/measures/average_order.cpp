#include "./average_order.hpp"

namespace inchworm::measures {

  average_order::average_order(params_t const &, config_t const &config, qmc_results_t &results)
     : config(config), average_order_ref(results.average_order) {
    average_order_ref = 0.0;
  }

  void average_order::accumulate(scalar_t) {
    average_order_ref += config.size();
    ++N;
  }

  void average_order::collect_results(mpi::communicator const &comm) {
    N = mpi::all_reduce(N, comm);

    // Reduce and normalize
    average_order_ref = mpi::all_reduce(average_order_ref, comm);
    average_order_ref = average_order_ref / N;
  }

} // namespace inchworm::measures
