#include "./U_frame.hpp"

namespace inchworm::measures {

  U_frame::U_frame(params_t const &params, qmc_data_t const &qmc_data_, container_set &results_)
     : qmc_data(qmc_data_), results(results_), average_sign_(0.0), U_frame_(qmc_data.u_frame) {
    U_frame_.reset();
  }

  void U_frame::accumulate(mc_weight_t sign) {
    U_frame_ += qmc_data.u_frame;
    average_sign_ += sign;
    //++count;
  }

  void U_frame::collect_results(mpi::communicator const &comm) {
    average_sign_ = mpi::all_reduce(average_sign_, comm);
    //U_frame_      = mpi::all_reduce(U_frame_, comm);  //???????

    //assign_frame_to_propagator(results.u_tau, U_frame_, qmc_data.inch_step);
  }

} // namespace inchworm::measures
