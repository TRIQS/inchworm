#include "./insert.hpp"

namespace inchworm::moves {

  scalar_t insert::try_move(config_t &config) {
    int n_fops = all_d_ops.size();
    auto d     = all_d_ops[rng(n_fops)];
    auto d_dag = all_d_dag_ops[rng(n_fops)];

    d.tau     = rng(params.tau_max);
    d_dag.tau = rng(params.tau_max);

    if (not config.try_insert(d_dag, d)) return 0.0;

    int N = config.size();
    return std::pow(params.tau_max * n_fops / N, 2);
  }

  scalar_t double_insert::try_move(config_t &config) {
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
