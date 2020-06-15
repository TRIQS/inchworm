#include "./u_frame.hpp"

namespace inchworm::measures {

  u_frame::u_frame(params_t const &, qmc_config_data_t const &qmc_config_data_, single_step_results_t &results_)
     : qmc_config_data(qmc_config_data_), results(results_) {}

  void u_frame::accumulate(scalar_t sign) {
    scalar_t s = sign / (qmc_config_data.w.loc);
    if (qmc_config_data.config.size() < MAX_ORDER) { // MAX_ORDER is just for printing purpose for now.
      auto const &[tbl0, m0] = qmc_config_data.u_partial[0];
      if (tbl0 != -1) results.u_expansion_order[qmc_config_data.config.size()] += s * m0(0, 0);
    }
    //if (qmc_config_data.config.size() == 2)
    //std::printf("%d ", qmc_config_data.config.size());
    for (int bl = 0; bl < results.u_frame.size(); bl++) {
      auto const &[tbl, m] = qmc_config_data.u_partial[bl];
      if (tbl != -1) results.u_frame[bl] += s * m;
    }

    // For normalizatoin purpose, we sample the zeroth order separatly:
    if (qmc_config_data.config.size() == 0)
      for (int bl = 0; bl < results.u_frame_0th_order.size(); bl++) {
        auto const &[tbl, m] = qmc_config_data.u_partial[bl];
        if (tbl != -1) results.u_frame_0th_order[bl] += s * m;
      }
  }

  void u_frame::collect_results(mpi::communicator const &comm) {
    results.u_frame_0th_order = mpi::all_reduce(results.u_frame_0th_order, comm);
    results.u_frame           = mpi::all_reduce(results.u_frame, comm);
    results.u_expansion_order = mpi::all_reduce(results.u_expansion_order, comm);
  }

} // namespace inchworm::measures
