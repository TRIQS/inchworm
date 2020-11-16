#include "./remove.hpp"

namespace inchworm::moves {

  scalar_t remove::try_move(config_t &config) {
    int N = config.size();
    if (N == 0) return 0;

    int idx     = rng(N);
    int idx_dag = rng(N);

    if (not config.try_erase(idx_dag, idx)) return 0;

    int n_fops = all_d_ops.size();
    return std::pow(N / (n_fops * params.tau_max), 2);
  }

  scalar_t double_remove::try_move(config_t &config) {
    int N = config.size();
    if (N == 0) return 0;

    int idx1     = rng(N);
    int idx1_dag = rng(N);
    int idx2     = rng(N);
    int idx2_dag = rng(N);

    if (not config.try_double_erase(idx1_dag, idx1, idx2_dag, idx2)) return 0;

    int n_fops = all_d_ops.size();
    return std::pow(N / (n_fops * params.tau_max), 4);
  }

} // namespace inchworm::moves
