#include "./u_frame.hpp"

namespace inchworm::measures {

  //u_frame::u_frame(params_t const &, qmc_config_data_t const &qmc_config_data_, container_set &results_)
  u_frame::u_frame(params_t const &, qmc_config_data_t const &qmc_config_data_, single_step_results_t &results_)
     : qmc_config_data(qmc_config_data_), results(results_) {}

  void u_frame::accumulate(scalar_t sign) {
    //int factor = 1;
    scalar_t s = sign / qmc_config_data.w.loc;
    //if (qmc_config_data.config.size() < MAX_ORDER) { results.u_expansion_order[qmc_config_data.config.size()] += 1 / qmc_config_data.w.loc; }
    //average_sign += s;
    if (qmc_config_data.config.size() == 0)
      for (int bl = 0; bl < results.zero_frame.size(); bl++) results.zero_frame[bl] += s * qmc_config_data.u_frame[bl];
    if (qmc_config_data.config.size() == 1)
      for (int bl = 0; bl < results.u_frame.size(); bl++) results.u_frame[bl] += s * qmc_config_data.u_frame[bl];
  }

  void u_frame::collect_results(mpi::communicator const &comm) {
    //results.zero_frame = mpi::all_reduce(results.zero_frame, comm);
    results.u_frame    = mpi::all_reduce(results.u_frame, comm);
    //results.u_expansion_order = mpi::all_reduce(results.u_expansion_order, comm);
    //for (auto &x : results.u_frame) x /= 10000;
  }

} // namespace inchworm::measures
