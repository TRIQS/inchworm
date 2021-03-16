#include "./frame.hpp"
#include "./../u_frame.hpp"

namespace inchworm::measures {

  frame::frame(params_t const &params, qmc_data_t const &qmc_data_, qmc_results_t &results)
     : verbosity(params.verbosity), qmc_data(qmc_data_), frame_ref(results.frame), frame_0th_order_ref(results.frame_0th_order) {}

  void frame::accumulate(scalar_t sign) {

    // We weight the Monte-Carlo by the frobenius norm of the current frame
    // This importance sampling factor has to be corrected in the measurement
    scalar_t s = sign / (qmc_data.weights.imp);

    for (int bl : range(frame_ref.size()))
      if (not qmc_data.frame[bl].empty()) frame_ref[bl] += s * qmc_data.frame[bl];

    // For normalization purpose, we sample the zeroth order separatly:
    if (qmc_data.config.size() == 0) frame_0th_order_ref += s * qmc_data.frame;

    // Perform an autocorrelation analysis on the trace
    acc << trace(qmc_data.frame);
  }

  void frame::collect_results(mpi::communicator const &comm) {
    frame_0th_order_ref = mpi::all_reduce(frame_0th_order_ref, comm);
    frame_ref           = mpi::all_reduce(frame_ref, comm);

    auto [errs, counts]  = acc.log_bin_errors_all_reduce(comm);

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
      std::printf("     auto_corr_time: %.3f\n", auto_corr_time);
    }
  }

} // namespace inchworm::measures
