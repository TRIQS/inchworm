#include "./autocorr.hpp"
#include "./../u_frame.hpp"

#include <h5/h5.hpp>

namespace inchworm::measures {

  autocorr::autocorr(params_t const &params, qmc_data_t const &qmc_data_, qmc_results_t &results_)
     : verbosity(params.verbosity), qmc_data(qmc_data_), results(results_) {}

  void autocorr::accumulate(scalar_t sign) {

    // We weight the Monte-Carlo by the frobenius norm of the current frame
    // This importance sampling factor has to be corrected in the measurement
    scalar_t s = sign / (qmc_data.weights.imp);

    //acc << s * trace(qmc_data.frame);
    acc << qmc_data.config.size();
  }

  void autocorr::collect_results(mpi::communicator const &comm) {

    auto errs = acc.log_bin_errors_mpi(comm);

    // Debug Prints
    if (comm.rank() == 0 and verbosity > 1) {
      std::cout << "errs avg_k: [";
      for (auto err : errs) std::cout << err << " ";
      std::cout << "]\n";
    }

    // Estimate auto-correlation time
    results.auto_corr_time = 0.0;
    if (comm.rank() == 0 && errs[0] > 0) results.auto_corr_time = std::max(0.0, tau_estimate_from_errors(errs[int(0.7 * errs.size())], errs[0]));
    mpi::broadcast(results.auto_corr_time, comm, 0);

    // Reset the accumulator
    acc = {0.0, -1};
  }

} // namespace inchworm::measures
