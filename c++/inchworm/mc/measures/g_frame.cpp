#include "./g_frame.hpp"

namespace inchworm::measures {

  g_frame::g_frame(params_t const &, qmc_config_data_t const &qmc_config_data_, single_step_results_t &results_)
     : qmc_config_data(qmc_config_data_), results(results_) {}

  void g_frame::accumulate(scalar_t sign) {

    scalar_t s = sign / (qmc_config_data.w.loc);

    size_t pert_order = qmc_config_data.config.size();
    if (pert_order >= results.expansion_order.size()) {
      size_t new_size = std::max(2 * results.expansion_order.size(), pert_order + 1);
      results.expansion_order.resize(new_size, 0);
      results.samples_expansion_order.resize(new_size, 0);
    }

    auto m0 = qmc_config_data.g_frame[0](0, 0);
    results.expansion_order[pert_order] += s * m0;
    results.samples_expansion_order[pert_order]++;
    results.measure_count++;
    results.average_k = results.average_k + pert_order;

    for (int bl = 0; bl < results.frame.size(); bl++) {
      auto m = qmc_config_data.g_frame[bl];
      results.frame[bl] += s * m;
    }

    // For normalizatoin purpose, we sample the zeroth order separatly:
    if (pert_order == 0) {
      for (int bl = 0; bl < results.frame_0th_order.size(); bl++) {
        auto m = qmc_config_data.g_frame[bl];
        results.frame_0th_order[bl] += s * m;
      }
    }
  }

  void g_frame::collect_results(mpi::communicator const &comm) {
    results.frame_0th_order = mpi::all_reduce(results.frame_0th_order, comm);
    results.frame           = mpi::all_reduce(results.frame, comm);
    results.measure_count   = mpi::all_reduce(results.measure_count, comm);
    results.average_k       = mpi::all_reduce(results.average_k, comm) / results.measure_count;

    auto max_size = mpi::all_reduce(results.expansion_order.size(), comm, MPI_MAX);
    results.expansion_order.resize(max_size, 0);
    results.expansion_order = mpi::all_reduce(results.expansion_order, comm);

    max_size = mpi::all_reduce(results.samples_expansion_order.size(), comm, MPI_MAX);
    results.samples_expansion_order.resize(max_size, 0.0);
    results.samples_expansion_order = mpi::all_reduce(results.samples_expansion_order, comm);
  }

} // namespace inchworm::measures
