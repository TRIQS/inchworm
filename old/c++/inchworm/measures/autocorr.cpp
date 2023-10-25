#include "./autocorr.hpp"

#include <h5/h5.hpp>

namespace inchworm::measures {

  autocorr::autocorr(params_t const &params, config_t const &config, frame_t const &frame, qmc_results_t &results)
     : verbosity(params.verbosity), config(config), curr_frame(frame), results(results) {}

  void autocorr::accumulate(scalar_t sign) {
    // We weight the Monte-Carlo by the frobenius norm of the current frame
    // This importance sampling factor has to be corrected in the measurement
    scalar_t s = sign / (config.imp_weight);

    log_acc << trace(s * curr_frame);
  }

  void autocorr::collect_results(mpi::communicator const &comm) {

    auto [errs, counts] = log_acc.log_bin_errors_all_reduce(comm);

    // Debug Prints
    if (comm.rank() == 0 and verbosity > 1) {
      std::cout << "errs norm(frame): [";
      for (auto err : errs) std::cout << err << " ";
      std::cout << "]\n";
    }

    // Estimate auto-correlation time
    results.auto_corr_time = 0.0;
    if (comm.rank() == 0 && errs[0] > 0) results.auto_corr_time = std::max(0.0, tau_estimate_from_errors(errs[int(0.7 * errs.size())], errs[0]));
    mpi::broadcast(results.auto_corr_time, comm, 0);

    // Reset the accumulator
    log_acc = {0.0, -1, 0};
  }

} // namespace inchworm::measures
