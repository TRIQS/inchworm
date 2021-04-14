#include "./insert.hpp"
#include "./gen_op_times.hpp"
#include "./../torus.hpp"

namespace inchworm::moves {

  scalar_t insert::try_config_update(config_t &config) {
    long n_bl = gf_struct.size();

    long bl      = rng(n_bl);
    long bl_size = gf_struct[bl].second;

    auto d     = solver.all_d_ops[bl][rng(bl_size)];
    auto d_dag = solver.all_d_dag_ops[bl][rng(bl_size)];

    auto [d_tau, d_dag_tau] = gen_op_times(rng, d, d_dag, config, params.tau_split, params.tau_max);
    d.tau                   = d_tau;
    d_dag.tau               = d_dag_tau;

    last_d     = d;
    last_d_dag = d_dag;

    double prop_prob = get_time_prop_prob(d, d_dag, config, params.tau_split, params.tau_max) / bl_size / bl_size / n_bl;

    if (not config.try_insert(d_dag, d)) return 0.0;

    double inv_prop_prob = 1.0 / config.size(bl) / config.size(bl) / n_bl;

    return inv_prop_prob / prop_prob;
  }

  scalar_t insert::accept() {

    // Gather tau statistic during warmup
    if (gather_tau_diff_stat) {
      for (auto const &op : {last_d, last_d_dag}) {
        double tau_diff = cyclic_difference(op, prop_config, params.tau_max);
        if (tau_diff < params.tau_max / 2.0) // left_width
          tau_diff_stat[op.bl](op.idx, op.dag, 0) << tau_diff;
        else // right_width
          tau_diff_stat[op.bl](op.idx, op.dag, 1) << params.tau_max - tau_diff;
      }
    }

    return base_move::accept();
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

    auto d1     = solver.all_d_ops[bl1][rng(bl1_size)];
    auto d1_dag = solver.all_d_dag_ops[bl1][rng(bl1_size)];
    auto d2     = solver.all_d_ops[bl2][rng(bl2_size)];
    auto d2_dag = solver.all_d_dag_ops[bl2][rng(bl2_size)];

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
