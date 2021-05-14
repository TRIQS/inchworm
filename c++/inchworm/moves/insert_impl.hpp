#pragma once

#include "insert.hpp"
#include "gen_op_times.hpp"
#include "./../torus.hpp"

namespace inchworm::moves {

  template <KIND Kind> scalar_t insert<Kind>::insert_single_op_pair(long bl, config_t &config) {

    long bl_size = gf_struct[bl].second;

    auto d     = solver.all_d_ops[bl][rng(bl_size)];
    auto d_dag = solver.all_d_dag_ops[bl][rng(bl_size)];

    auto [d_tau, d_dag_tau] = gen_op_times(rng, d, d_dag, config, params.tau_split, params.tau_max);
    d.tau                   = d_tau;
    d_dag.tau               = d_dag_tau;

    last_d     = d;
    last_d_dag = d_dag;

    double prop_prob = get_time_prop_prob(d, d_dag, config, params.tau_split, params.tau_max) / bl_size / bl_size;

    if (not config.try_insert(d_dag, d)) return 0.0;

    double inv_prop_prob = 1.0 / config.size(bl) / config.size(bl);

    return inv_prop_prob / prop_prob;
  }

  template <KIND Kind> scalar_t insert<Kind>::insert_double_op_pair(long bl1, long bl2, config_t &config) {

    long bl1_size = gf_struct[bl1].second;
    long bl2_size = gf_struct[bl2].second;

    auto d1     = solver.all_d_ops[bl1][rng(bl1_size)];
    auto d1_dag = solver.all_d_dag_ops[bl1][rng(bl1_size)];
    auto d2     = solver.all_d_ops[bl2][rng(bl2_size)];
    auto d2_dag = solver.all_d_dag_ops[bl2][rng(bl2_size)];

    d1.tau                   = rng(params.tau_max);
    d1_dag.tau               = rng(params.tau_max);
    d2.tau                   = rng(params.tau_max);
    d2_dag.tau               = rng(params.tau_max);

    last_d     = d2;
    last_d_dag = d2_dag;

    double prop_prob = std::pow(1.0 / params.tau_max, 4) / bl1_size / bl1_size / bl2_size / bl2_size;

    if (not config.try_double_insert(d1_dag, d1, d2_dag, d2)) return 0.0;

    double inv_prop_prob = 1.0 / config.size(bl1) / config.size(bl1) / config.size(bl2) / config.size(bl2);

    return inv_prop_prob / prop_prob;
  }

  template <KIND Kind> scalar_t insert<Kind>::try_config_update(config_t &config) {

    long n_bl = gf_struct.size();

    if constexpr (Kind == KIND::Single) {

      return insert_single_op_pair(rng(n_bl), config);

    } else { // Double Insert

      long bl1 = rng(n_bl);
      long bl2 = (Kind == KIND::DoubleEqBl) ? bl1 : rng(n_bl);

      //return insert_double_op_pair(bl1, bl2, config);

      double t_ratio1 = insert_single_op_pair(bl1, config);
      if (t_ratio1 == 0.0) return 0.0;

      double t_ratio2 = insert_single_op_pair(bl2, config);

      return t_ratio1 * t_ratio2;
    }
  }

  template <KIND Kind> scalar_t insert<Kind>::accept() {

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

} // namespace inchworm::moves
