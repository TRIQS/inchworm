#pragma once

#include "./remove.hpp"
#include "./gen_op_times.hpp"

namespace inchworm::moves {

  template <KIND Kind> scalar_t remove<Kind>::remove_single_op_pair(long bl, config_t &config) {

    if (config.size(bl) == 0) return 0;

    long idx     = rng(config.size(bl));
    long idx_dag = rng(config.size(bl));

    fop_t d     = config.d_bl_list[bl][idx];
    fop_t d_dag = config.d_dag_bl_list[bl][idx];

    double prop_prob = 1.0 / config.size(bl) / config.size(bl);

    // Never accept removal from trivial configurations
    if (params.tau_split > 0.0 && config.size() == 1)
      if ((d.tau < params.tau_split && d_dag.tau < params.tau_split) ||
          (d.tau > params.tau_split && d_dag.tau > params.tau_split)) return 0;

    if (not config.try_erase(bl, idx_dag, idx)) return 0;

    long bl_size         = gf_struct[bl].second;
    double inv_prop_prob = get_time_prop_prob(d, d_dag, config, params.tau_split, params.tau_max) / bl_size / bl_size;

    return inv_prop_prob / prop_prob;
  }

  template <KIND Kind> scalar_t remove<Kind>::remove_double_op_pair(long bl1, long bl2, config_t &config) {

    if (config.size(bl1) == 0) return 0;
    if (config.size(bl2) == 0) return 0;

    long idx1     = rng(config.size(bl1));
    long idx1_dag = rng(config.size(bl1));
    long idx2     = rng(config.size(bl2));
    long idx2_dag = rng(config.size(bl2));

    fop_t d1     = config.d_bl_list[bl1][idx1];
    fop_t d1_dag = config.d_dag_bl_list[bl1][idx1];
    fop_t d2     = config.d_bl_list[bl2][idx2];
    fop_t d2_dag = config.d_dag_bl_list[bl2][idx2];

    double prop_prob = 1.0 / config.size(bl1) / config.size(bl1) / config.size(bl2) / config.size(bl2);

    if (not config.try_double_erase(bl1, idx1_dag, idx1, bl2, idx2_dag, idx2)) return 0;

    long bl1_size         = gf_struct[bl1].second;
    long bl2_size         = gf_struct[bl2].second;
    double inv_prop_prob  = std::pow(1.0 / params.tau_max, 4) / bl1_size / bl1_size / bl2_size / bl2_size;

    return inv_prop_prob / prop_prob;
  }

  template <KIND Kind> scalar_t remove<Kind>::try_config_update(config_t &config) {

    long n_bl = gf_struct.size();

    if constexpr (Kind == KIND::Single) {

      return remove_single_op_pair(rng(n_bl), config);

    } else { // Double Insert

      long bl1 = rng(n_bl);
      long bl2 = (Kind == KIND::DoubleEqBl) ? bl1 : rng(n_bl);

      if (bl1 == bl2 and config.size(bl1) < 2) return 0;

      //return remove_double_op_pair(bl1, bl2, config);

      double t_ratio1 = remove_single_op_pair(bl1, config);
      if (t_ratio1 == 0.0) return 0.0;

      double t_ratio2 = remove_single_op_pair(bl2, config);

      return t_ratio1 * t_ratio2;
    }
  }

} // namespace inchworm::moves
