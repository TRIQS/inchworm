#include "./g_frame.hpp"

namespace inchworm::measures {

  g_frame::g_frame(params_t const &, qmc_config_data_t const &qmc_config_data_, single_step_results_t &results_)
     : qmc_config_data(qmc_config_data_), results(results_) {}

  void g_frame::accumulate(scalar_t sign) {
    scalar_t s = sign / (qmc_config_data.w.loc);
    if (qmc_config_data.config.size() < MAX_ORDER) { // MAX_ORDER is just for printing purpose for now.
      auto m0 = qmc_config_data.g_frame[0](0,0);
      results.expansion_order[qmc_config_data.config.size()] += s * m0;
    }
    //if (qmc_config_data.config.size() == 2)
    //std::printf("%d ", qmc_config_data.config.size());
    for (int bl = 0; bl < results.frame.size(); bl++) {
      auto m = qmc_config_data.g_frame[bl];
      results.frame[bl] += s * m;
    }

    // For normalizatoin purpose, we sample the zeroth order separatly:
    if (qmc_config_data.config.size() == 0)
      for (int bl = 0; bl < results.frame_0th_order.size(); bl++) {
        auto m = qmc_config_data.g_frame[bl];
        results.frame_0th_order[bl] += s * m;
      }
  }

  void g_frame::collect_results(mpi::communicator const &comm) {
    results.frame_0th_order = mpi::all_reduce(results.frame_0th_order, comm);
    results.frame           = mpi::all_reduce(results.frame, comm);
    results.expansion_order = mpi::all_reduce(results.expansion_order, comm);
  }

} // namespace inchworm::measures
