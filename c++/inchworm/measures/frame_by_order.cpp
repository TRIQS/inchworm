#include "./frame_by_order.hpp"

namespace inchworm::measures {

  frame_by_order::frame_by_order(params_t const &, qmc_data_t const &qmc_data_, qmc_results_t &results_)
     : qmc_data(qmc_data_), frame_by_order_ref(results_.frame_by_order), zero_frame(results_.frame) {
    for (auto &bl : zero_frame) bl = 0.;
    frame_by_order_ref = std::vector<frame_t>(10, zero_frame);
  }

  void frame_by_order::accumulate(scalar_t sign) {

    // Resize the vector as necessary
    size_t pert_order = qmc_data.config.size();
    if (pert_order >= frame_by_order_ref->size()) {
      size_t new_size = std::max(2 * frame_by_order_ref->size(), pert_order + 1);
      frame_by_order_ref->resize(new_size, zero_frame);
    }

    scalar_t s = sign / (qmc_data.weights.imp);
    for (int bl : range(qmc_data.frame.size())) {
      if (not qmc_data.frame[bl].empty()) { (*frame_by_order_ref)[pert_order][bl] += s * qmc_data.frame[bl]; }
    }
  }

  void frame_by_order::collect_results(mpi::communicator const &comm) {

    // Make sure that all mpi threads have an equally sized histogram
    auto max_size = mpi::all_reduce(frame_by_order_ref->size(), comm, MPI_MAX);
    frame_by_order_ref->resize(max_size, zero_frame);

    // Reduce histogram over all mpi threads
    frame_by_order_ref = mpi::all_reduce(*frame_by_order_ref, comm);
  }

} // namespace inchworm::measures
