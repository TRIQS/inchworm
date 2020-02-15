#include "./sign.hpp"

namespace inchworm::measures {

  sign::sign(params_t const &, qmc_config_data_t const &qmc_config_data_, single_step_results_t &results_)
     : qmc_config_data(qmc_config_data_), results(results_) {}

  void sign::accumulate(scalar_t sign) {
    double atomic_reweighting = trace(qmc_config_data.u_frame) / qmc_config_data.w.loc;
    average_sign += sign * atomic_reweighting;
    z += std::abs(atomic_reweighting);
  }

  void sign::collect_results(mpi::communicator const &comm) {
    average_sign         = mpi::all_reduce(average_sign, comm);
    z                    = mpi::all_reduce(z, comm);
    results.average_sign = average_sign / z;
    std::printf("z=%d average_sign=%d\n", z, average_sign);
  }

} // namespace inchworm::measures
