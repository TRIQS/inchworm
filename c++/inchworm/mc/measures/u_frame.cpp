#include "./u_frame.hpp"

namespace inchworm::measures {

  //u_frame::u_frame(params_t const &, qmc_config_data_t const &qmc_config_data_, container_set &results_)
  u_frame::u_frame(params_t const &, qmc_config_data_t const &qmc_config_data_, single_step_results_t &results_)
     : qmc_config_data(qmc_config_data_), results(results_) {}

  void u_frame::accumulate(scalar_t sign) {
    //if (qmc_config_data.config.size() == 1) {
      for (int bl = 0; bl < results.u_frame.size(); bl++) results.u_frame[bl] += qmc_config_data.u_frame[bl];
      average_sign += sign;
    //}
  }

  void u_frame::collect_results(mpi::communicator const &comm) {
    average_sign    = mpi::all_reduce(average_sign, comm);
    results.u_frame = mpi::all_reduce(results.u_frame, comm);

    //std::printf("average sign: % 3.4f \n", average_sign);
    for (auto &x : results.u_frame) x /= average_sign;
    //assign_frame_to_propagator(results.u_tau, u_frame, qmc_config_data.inch_step);
  }

} // namespace inchworm::measures
