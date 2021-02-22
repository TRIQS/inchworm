#include "./remove.hpp"

namespace inchworm::moves {

  scalar_t remove::try_config_update(config_t &config) {
    int OldSize = config.size();
    if (OldSize == 0) return 0;

    int idx     = rng(OldSize);
    int idx_dag = rng(OldSize);

    if (not config.try_erase(idx_dag, idx)) return 0;

    int n_fops = all_d_ops.size();

    // Account for special treatment of insertion into empty configuration
    if (OldSize == 1 && (params.tau_split != 0.0)) {
      auto dtau = params.tau_max - params.tau_split;
      return 1.0 * OldSize * OldSize / (2.0 * params.tau_split * dtau * n_fops * n_fops);
    } else {
      return std::pow(OldSize / (n_fops * params.tau_max), 2);
    }
  }

  scalar_t double_remove::try_config_update(config_t &config) {
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
