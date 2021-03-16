#include "./insert.hpp"

namespace inchworm::moves {

  scalar_t insert::try_config_update(config_t &config) {
    long n_bl = gf_struct.size();

    long bl      = rng(n_bl);
    long bl_size = gf_struct[bl].second;

    auto d     = all_d_ops[bl][rng(bl_size)];
    auto d_dag = all_d_dag_ops[bl][rng(bl_size)];

    scalar_t t_ratio;
    if (config.size() == 0 && (params.tau_split != 0.0)) {
      // Special case of empty configuration for inchworm
      // Make sure that we choose tau values on seperate
      // sides of the split point

      auto dtau = params.tau_max - params.tau_split;
      d.tau     = rng(params.tau_split);
      d_dag.tau = params.tau_split + rng(dtau);
      if (rng(2)) std::swap(d.tau, d_dag.tau);
      t_ratio = 2.0 * params.tau_split * dtau * bl_size * bl_size;
    } else {
      d.tau           = rng(params.tau_max);
      d_dag.tau       = rng(params.tau_max);
      long new_nop_bl = config.size(bl) + 1;
      t_ratio         = std::pow(params.tau_max * bl_size / new_nop_bl, 2);
    }

    if (not config.try_insert(d_dag, d)) return 0.0;

    return t_ratio;
  }

  scalar_t double_insert::try_config_update(config_t &config) {
    long n_bl = gf_struct.size();

    long bl1, bl2;
    if (equal_blocks) {
      bl1 = rng(n_bl);
      bl2 = bl1;
    } else {
      bl1 = rng(n_bl);
      bl2 = rng(n_bl);
    }

    long bl1_size = gf_struct[bl1].second;
    long bl2_size = gf_struct[bl2].second;

    auto d1     = all_d_ops[bl1][rng(bl1_size)];
    auto d1_dag = all_d_dag_ops[bl1][rng(bl1_size)];
    auto d2     = all_d_ops[bl2][rng(bl2_size)];
    auto d2_dag = all_d_dag_ops[bl2][rng(bl2_size)];

    d1.tau     = rng(params.tau_max);
    d1_dag.tau = rng(params.tau_max);
    d2.tau     = rng(params.tau_max);
    d2_dag.tau = rng(params.tau_max);

    if (not config.try_double_insert(d1_dag, d1, d2_dag, d2)) return 0;

    long new_nop_bl1 = config.d_bl_list[bl1].size();
    long new_nop_bl2 = config.d_bl_list[bl2].size();
    return std::pow(params.tau_max * params.tau_max * bl1_size * bl2_size / (new_nop_bl1 * new_nop_bl2), 2);
  }

} // namespace inchworm::moves
