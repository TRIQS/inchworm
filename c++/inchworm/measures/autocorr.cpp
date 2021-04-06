#include "./autocorr.hpp"

#include <h5/h5.hpp>

namespace inchworm::measures {

  autocorr::autocorr(params_t const &params, config_t const &config, qmc_results_t &results)
     : verbosity(params.verbosity), config(config), results(results) {}

  void autocorr::accumulate(scalar_t sign) { acc << config.size(); }

  void autocorr::collect_results(mpi::communicator const &comm) {

    auto [errs, counts] = acc.log_bin_errors_all_reduce(comm);

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
