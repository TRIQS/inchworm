#include "./u_frame.hpp"

namespace inchworm::measures {

  u_frame::u_frame(params_t const &params, qmc_config_data_t const &qmc_config_data_, container_set &results_)
     : qmc_config_data(qmc_config_data_), results(results_), average_sign_(0.0), u_frame_(qmc_config_data.u_frame) {
    //u_frame.reset();
  }

  void u_frame::accumulate(scalar_t sign) {
    //u_frame += qmc_config_data.u_frame;
    for (int bl = 0; bl < u_frame_.size(); bl++) u_frame_[bl] += qmc_config_data.u_frame[bl];
    average_sign_ += sign;
    //++count;
  }

  void u_frame::collect_results(mpi::communicator const &comm) {
    average_sign_ = mpi::all_reduce(average_sign_, comm);
    u_frame_      = mpi::all_reduce(u_frame_, comm);  //???????

    //assign_frame_to_propagator(results.u_tau, u_frame, qmc_config_data.inch_step);
  }

} // namespace inchworm::measures
