#include "./average_order.hpp"

namespace inchworm::measures {

  average_order::average_order(params_t const &, qmc_data_t const &qmc_data_, qmc_results_t &results)
     : qmc_data(qmc_data_), average_order_ref(results.average_order) {
    average_order_ref = 0.0;
  }

  void average_order::accumulate(scalar_t) {
    average_order_ref += qmc_data.config.size();
    ++N;
  }

  void average_order::collect_results(mpi::communicator const &comm) {
    N = mpi::all_reduce(N, comm);

    // Reduce and normalize
    average_order_ref = mpi::all_reduce(average_order_ref, comm);
    average_order_ref = average_order_ref / N;
  }

} // namespace inchworm::measures
