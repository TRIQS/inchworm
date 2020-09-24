#include "./u_frame.hpp"

namespace inchworm::measures {

  u_frame::u_frame(params_t const &, qmc_config_data_t const &qmc_config_data_, single_step_results_t &results_)
     : qmc_config_data(qmc_config_data_), results(results_) {}

  void u_frame::accumulate(scalar_t sign) {

    scalar_t s = sign / (qmc_config_data.w.loc);

    size_t pert_order = qmc_config_data.config.size();
    if (pert_order >= results.expansion_order.size()) {
      size_t new_size = std::max(2 * results.expansion_order.size(), pert_order + 1);
      results.expansion_order.resize(new_size, 0);
      results.samples_expansion_order.resize(new_size, 0);
    }

    auto const &[tbl0, m0] = qmc_config_data.u_partial[0];
    if (tbl0 != -1) results.expansion_order[pert_order] += s * m0(0, 0);
    results.samples_expansion_order[pert_order]++;
    results.measure_count++;
    results.average_k += pert_order;

    for (int bl = 0; bl < results.frame.size(); bl++) {
      auto const &[tbl, m] = qmc_config_data.u_partial[bl];
      if (tbl != -1) results.frame[bl] += s * m;
    }

    // For normalizatoin purpose, we sample the zeroth order separatly:
    if (pert_order == 0) {
      for (int bl = 0; bl < results.frame_0th_order.size(); bl++) {
        auto const &[tbl, m] = qmc_config_data.u_partial[bl];
        if (tbl != -1) results.frame_0th_order[bl] += s * m;
      }
    }
  }

  void u_frame::collect_results(mpi::communicator const &comm) {
    results.frame_0th_order         = mpi::all_reduce(results.frame_0th_order, comm);
    results.frame                   = mpi::all_reduce(results.frame, comm);
    results.measure_count           = mpi::all_reduce(results.measure_count, comm);
    results.average_k               = mpi::all_reduce(results.average_k, comm) / results.measure_count;
    results.expansion_order         = mpi::all_reduce(results.expansion_order, comm);
    results.samples_expansion_order = mpi::all_reduce(results.samples_expansion_order, comm);
  }

} // namespace inchworm::measures
