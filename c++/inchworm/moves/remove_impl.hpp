// Copyright (c) 2021--present, The Simons Foundation
// This file is part of inchworm and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once

#include "./remove.hpp"
#include "./gen_op_times.hpp"

namespace inchworm::moves {

  template <KIND Kind> double remove<Kind>::remove_single_op_pair(long bl, config_t &config) {

    if (config.size(bl) == 0) return 0;

    long idx     = rng(config.size(bl));
    long idx_dag = rng(config.size(bl));

    fop_t d     = config.d_bl_list[bl][idx];
    fop_t d_dag = config.d_dag_bl_list[bl][idx_dag];

    double prop_prob = 1.0 / config.size(bl) / config.size(bl);

    // Never accept removal from trivial configurations of first order
    // Note: Special insertion rule will make sure that 1st order is non-trivial
    if (params.tau_split > 0.0 && config.size() == 1)
      if ((d.tau < params.tau_split && d_dag.tau < params.tau_split) || (d.tau > params.tau_split && d_dag.tau > params.tau_split)) return 0;

    if (not config.try_erase(bl, idx_dag, idx)) return 0;

    long bl_size         = gf_struct[bl].second;
    double inv_prop_prob = get_time_prop_prob(d, d_dag, config, params.tau_split, params.tau_max) / bl_size / bl_size;

    return inv_prop_prob / prop_prob;
  }

  template <KIND Kind> double remove<Kind>::try_config_update(config_t &config) {

    long n_bl = gf_struct.size();

    if constexpr (Kind == KIND::Single) {

      return remove_single_op_pair(rng(n_bl), config);

    } else { // Double Insert

      long bl1 = rng(n_bl);
      long bl2 = (Kind == KIND::DoubleEqBl) ? bl1 : rng(n_bl);

      if (bl1 == bl2 and config.size(bl1) < 2) return 0;

      double t_ratio1 = remove_single_op_pair(bl1, config);
      if (t_ratio1 == 0.0) return 0.0;

      double t_ratio2 = remove_single_op_pair(bl2, config);

      return t_ratio1 * t_ratio2;
    }
  }

} // namespace inchworm::moves
