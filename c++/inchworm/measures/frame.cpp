#include "./frame.hpp"
#include "./../u_frame.hpp"

namespace inchworm::measures {

  frame::frame(params_t const &params, config_t const &config, frame_t const &frame, qmc_results_t &results)
     : verbosity(params.verbosity),
       config(config),
       curr_frame(frame),
       acc_frame(results.frame),
       acc_frame_zeroth_order(results.frame_zeroth_order),
       errs_frame(results.errs_frame),
       lin_acc(frame.size(), accumulator<scalar_t>{0.0, 0, 1000}) {
    for (auto &bl : acc_frame) bl = 0.;
    for (auto &bl : acc_frame_zeroth_order) bl = 0.;
  }

  void frame::accumulate(scalar_t sign) {

    // We weight the Monte-Carlo by the frobenius norm of the current frame
    // This importance sampling factor has to be corrected in the measurement
    scalar_t s = sign / (config.imp_weight);

    // Perform sampling and error analysis
    for (int bl : range(curr_frame.size()))
      if (not curr_frame[bl].empty()) {
        acc_frame[bl] += s * curr_frame[bl];
        lin_acc[bl] << s * curr_frame[bl](0, 0);
      } else {
        lin_acc[bl] << 0.0;
      }

    // For normalization purpose, we sample the zeroth order separatly:
    if (config.size() == 0) acc_frame_zeroth_order += s * curr_frame;

    // Perform an autocorrelation analysis on the trace
    log_acc << trace(s * curr_frame);

    ++N_samples;
  }

  void frame::collect_results(mpi::communicator const &comm) {
    acc_frame_zeroth_order = mpi::all_reduce(acc_frame_zeroth_order, comm);
    acc_frame              = mpi::all_reduce(acc_frame, comm);
    N_samples              = mpi::all_reduce(N_samples, comm);

    // Estimate error of the frame[0](0,0) component
    for (int bl : range(curr_frame.size())) errs_frame[bl] = N_samples * std::get<1>(mean_and_err_mpi(comm, lin_acc[bl].linear_bins()));

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
    lin_acc = {size_t(curr_frame.size()), accumulator<scalar_t>{0.0, 0, 1000}};
  }

} // namespace inchworm::measures
