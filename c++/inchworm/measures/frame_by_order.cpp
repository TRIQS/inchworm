#include "./frame_by_order.hpp"

namespace inchworm::measures {

  frame_by_order::frame_by_order(params_t const &, config_t const &config, frame_t const &frame, qmc_results_t &results)
     : config(config),
       frame_(frame),
       acc_frame_by_order(results.frame_by_order),
       err_frame_by_order(results.err_frame_by_order),
       zero_frame(results.frame) {
    for (auto &bl : zero_frame) bl = 0.;
  }

  void frame_by_order::accumulate(scalar_t sign) {

    // Resize the vector as necessary
    size_t pert_order = config.size();
    if (pert_order >= acc_frame_by_order.size()) {
      size_t new_size = std::max(2 * acc_frame_by_order.size(), pert_order + 1);
      acc_frame_by_order.resize(new_size, zero_frame);
      lin_acc_by_order.resize(new_size, empty_lin_acc);
    }

    scalar_t s = sign / (config.imp_weight);
    for (int bl : range(frame_.size())) {
      if (not frame_[bl].empty()) { acc_frame_by_order[pert_order][bl] += s * frame_[bl]; }
    }
    if (not frame_[0].empty()) lin_acc_by_order[pert_order] << s * frame_[0](0, 0);
  }

  void frame_by_order::collect_results(mpi::communicator const &comm) {

    // Make sure that all mpi threads have an equally sized histogram
    auto max_size = mpi::all_reduce(acc_frame_by_order.size(), comm, MPI_MAX);
    acc_frame_by_order.resize(max_size, zero_frame);
    lin_acc_by_order.resize(max_size, empty_lin_acc);

    // Reduce histogram over all mpi threads
    acc_frame_by_order = mpi::all_reduce(acc_frame_by_order, comm);

    // Estimate the errors of all the frame[k][0](0,0) components
    err_frame_by_order.clear();
    for (auto &lin_acc : lin_acc_by_order) {
      err_frame_by_order.push_back(std::get<1>(mean_and_err_mpi(comm, lin_acc.linear_bins())));
      lin_acc = empty_lin_acc;
    }
  }

} // namespace inchworm::measures
