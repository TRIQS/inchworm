#include "./sign.hpp"

namespace inchworm::measures {

  sign::sign(params_t const &params, qmc_config_t const &qmc_config_, container_set &results_)
     : qmc_config(qmc_config_), results(results_), average_sign_(0.0), count(0) {}

  void sign::accumulate(mc_weight_t sign) {
    average_sign_ += sign;
    ++count;
  }

  void sign::collect_results(mpi::communicator const &comm) {
    average_sign_ = mpi::all_reduce(average_sign_, comm);
    count         = mpi::all_reduce(count, comm);
    average_sign_ = average_sign_ / count;

    results.average_sign[qmc_config.inch_step] = average_sign_;
  }

} // namespace inchworm::measures
