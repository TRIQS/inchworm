#include "./insert.hpp"

namespace inchworm::moves {

  scalar_t insert::try_config_update(config_t &config) {
    int n_fops = all_d_ops.size();
    auto d     = all_d_ops[rng(n_fops)];
    auto d_dag = all_d_dag_ops[rng(n_fops)];

    scalar_t t_ratio;
    int NewSize = config.size() + 1;
    if (NewSize == 1 && (params.tau_split != 0.0)) {
      // Special case of empty configuration for inchworm
      // Make sure that we choose tau values on seperate
      // sides of the split point

      auto dtau = params.tau_max - params.tau_split;
      d.tau     = rng(params.tau_split);
      d_dag.tau = params.tau_split + rng(dtau);
      if (rng(2)) std::swap(d.tau, d_dag.tau);
      t_ratio = 2.0 * params.tau_split * dtau * n_fops * n_fops / NewSize / NewSize;
    } else {
      d.tau     = rng(params.tau_max);
      d_dag.tau = rng(params.tau_max);
      t_ratio   = std::pow(params.tau_max * n_fops / NewSize, 2);
    }

    if (not config.try_insert(d_dag, d)) return 0.0;

    return t_ratio;
  }

  scalar_t double_insert::try_config_update(config_t &config) {
    int n_fops  = all_d_ops.size();
    auto d1     = all_d_ops[rng(n_fops)];
    auto d1_dag = all_d_dag_ops[rng(n_fops)];
    auto d2     = all_d_ops[rng(n_fops)];
    auto d2_dag = all_d_dag_ops[rng(n_fops)];

    d1.tau     = rng(params.tau_max);
    d1_dag.tau = rng(params.tau_max);
    d2.tau     = rng(params.tau_max);
    d2_dag.tau = rng(params.tau_max);

    if (not config.try_double_insert(d1_dag, d1, d2_dag, d2)) return 0;

    int N = config.size();
    return std::pow(params.tau_max * n_fops / N, 4);
  }

} // namespace inchworm::moves
