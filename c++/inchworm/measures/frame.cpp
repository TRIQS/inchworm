#include "./frame.hpp"

namespace inchworm::measures {

  frame::frame(params_t const &, qmc_data_t const &qmc_data_, qmc_results_t &results)
     : qmc_data(qmc_data_), frame_ref(results.frame), frame_0th_order_ref(results.frame_0th_order) {}

  void frame::accumulate(scalar_t sign) {

    // We weight the Monte-Carlo by the frobenius norm of the current frame
    // This importance sampling factor has to be corrected in the measurement
    scalar_t s = sign / (qmc_data.weights.imp);

    for (int bl : range(frame_ref.size()))
      if (not qmc_data.frame[bl].empty()) frame_ref[bl] += s * qmc_data.frame[bl];

    // For normalization purpose, we sample the zeroth order separatly:
    if (qmc_data.config.size() == 0) frame_0th_order_ref += s * qmc_data.frame;
  }

  void frame::collect_results(mpi::communicator const &comm) {
    frame_0th_order_ref = mpi::all_reduce(frame_0th_order_ref, comm);
    frame_ref           = mpi::all_reduce(frame_ref, comm);
  }

} // namespace inchworm::measures
