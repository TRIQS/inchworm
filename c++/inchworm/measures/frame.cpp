#include "./frame.hpp"
#include "./../u_frame.hpp"

namespace inchworm::measures {

  frame::frame(params_t const &params, config_t const &config, frame_t const &frame, qmc_results_t &results)
     : verbosity(params.verbosity),
       config(config),
       frame_(frame),
       acc_frame(results.frame),
       err_frame(results.err_frame),
       acc_frame_0th_order(results.frame_0th_order) {}

  void frame::accumulate(scalar_t sign) {

    // We weight the Monte-Carlo by the frobenius norm of the current frame
    // This importance sampling factor has to be corrected in the measurement
    scalar_t s = sign / (config.imp_weight);

    for (int bl : range(frame_.size()))
      if (not frame_[bl].empty()) acc_frame[bl] += s * frame_[bl];

    // For normalization purpose, we sample the zeroth order separatly:
    if (config.size() == 0) acc_frame_0th_order += s * frame_;

    // Perform an autocorrelation analysis on the trace
    log_acc << trace(s * frame_);

    // Perform an error analysis on the [0](0,0) component
    if (not frame_[0].empty()) lin_acc << s * frame_[0](0, 0);
  }

  void frame::collect_results(mpi::communicator const &comm) {
    acc_frame_0th_order = mpi::all_reduce(acc_frame_0th_order, comm);
    acc_frame           = mpi::all_reduce(acc_frame, comm);

    // Estimate error of the frame[0](0,0) component
    err_frame = std::get<1>(mean_and_err_mpi(comm, lin_acc.linear_bins()));

    auto [errs, counts] = log_acc.log_bin_errors_all_reduce(comm);

    // Debug Prints
    if (comm.rank() == 0 and verbosity > 1) {
      std::cout << "errs Tr(frame): [";
      for (auto err : errs) std::cout << err << " ";
      std::cout << "]\n";
    }

    // Estimate auto-correlation time
    if (comm.rank() == 0) {
      double auto_corr_time = 0.0;
      if (errs[0] > 0) auto_corr_time = std::max(0.0, tau_estimate_from_errors(errs[int(0.7 * errs.size())], errs[0]));
      std::printf("     autocorr: %.3f\n", auto_corr_time);
    }

    // Reset the accumulators
    log_acc = {0.0, -1, 0};
    lin_acc = {0.0, 0, 1000};
  }

} // namespace inchworm::measures
