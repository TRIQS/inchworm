#include "./order_histogram.hpp"

namespace inchworm::measures {

  order_histogram::order_histogram(params_t const &params, config_t const &config, qmc_results_t &results)
     : config(config), order_histogram_ref(results.order_histogram) {
    order_histogram_ref = std::vector<double>(1 + params.max_order.value_or(10), 0.0);
  }

  void order_histogram::accumulate(scalar_t) {
    int k = config.size();
    while (k >= order_histogram_ref.size()) order_histogram_ref.resize(2 * order_histogram_ref.size());
    order_histogram_ref[k] += 1.;
    ++N;
  }

  void order_histogram::collect_results(mpi::communicator const &comm) {
    N = mpi::all_reduce(N, comm);

    // Make sure that all mpi threads have an equally sized histogram
    auto max_size = mpi::all_reduce(order_histogram_ref.size(), comm, MPI_MAX);
    order_histogram_ref.resize(max_size, 0.0);

    // Reduce histogram over all mpi threads
    order_histogram_ref = mpi::all_reduce(order_histogram_ref, comm);
    for (auto &h_k : order_histogram_ref) h_k = h_k / N;
  }

} // namespace inchworm::measures
