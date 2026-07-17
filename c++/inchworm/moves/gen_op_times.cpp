// Copyright (c) 2021--present, The Simons Foundation
// This file is part of inchworm and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "gen_op_times.hpp"

#include "./../distributions.hpp"

namespace inchworm::moves {

  std::array<double, 2> gen_op_times(triqs::mc_tools::random_generator &rng, fop_t const &d, fop_t const &d_dag, config_t const &config,
                                     double tau_split, double tau_max) {
    EXPECTS(d.bl == d_dag.bl);

    if (tau_split == 0.0) { // ====== CTHyb sampling with bare propagator ======

      return {rng(tau_max), rng(tau_max)};

    } else { // ============ Inchworm Sampling ===============

      if (config.size() == 0) { // ------- Empty Config -------

        // Make sure that we choose tau values on seperate sides of the split points at zero and tau_split
        double dtau = tau_max - tau_split;

        if (rng(2))
          return {double_icdf(rng(), d.left_width, d.right_width, tau_split),
                  tau_split + double_icdf(rng(), d_dag.left_width, d_dag.right_width, dtau)};
        else
          return {tau_split + double_icdf(rng(), d.left_width, d.right_width, dtau),
                  double_icdf(rng(), d_dag.left_width, d_dag.right_width, tau_split)};

      } else { // ------ Finite size config ------

        return {get_close_time(rng, d, config.split_times, tau_max), get_close_time(rng, d_dag, config.split_times, tau_max)};
      }
    }
  }

  double get_time_prop_prob(fop_t const &d, fop_t const &d_dag, config_t const &config, double tau_split, double tau_max) {
    EXPECTS(d.bl == d_dag.bl);

    if (tau_split == 0.0) { // ----- CTHyb sampling with bare propagator

      return 1.0 / tau_max / tau_max;

    } else { // ----- Inchworm Sampling

      if (config.size() == 0) { // Account for special treatment of empty config

        // Account for insertion around tau_split and zero
        double dtau = tau_max - tau_split;

        if (d.tau < tau_split)
          return double_pdf(d.tau, d.left_width, d.right_width, tau_split)
             * double_pdf(d_dag.tau - tau_split, d_dag.left_width, d_dag.right_width, dtau) / 2.0;
        else
          return double_pdf(d.tau - tau_split, d.left_width, d.right_width, dtau)
             * double_pdf(d_dag.tau, d_dag.left_width, d_dag.right_width, tau_split) / 2.0;

      } else { // ------ Finite block size ------

        return get_prob(d, config.split_times, tau_max) * get_prob(d_dag, config.split_times, tau_max);
      }
    }
  }

} // namespace inchworm::moves
