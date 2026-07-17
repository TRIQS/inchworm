// Copyright (c) 2020--present, The Simons Foundation
// This file is part of inchworm and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once

#include "base_move.hpp"

#include <numeric>

namespace inchworm::moves {

  template <KIND Kind> struct insert : public base_move {

    using acc_t = std::vector<nda::array<lin_binning<double>, 3>>;

    insert(config_t &config, frame_t &frame, qmc_params_t const &params, solver_core const &solver, triqs::mc_tools::random_generator &rng,
           acc_t &tau_diff_stat)
       : base_move(config, frame, params, solver, rng), tau_diff_stat(tau_diff_stat) {}

    inline std::string name() const override {
      switch (Kind) {
        case KIND::Single: return "Insert";
        case KIND::Double: return "Double Insert";
        case KIND::DoubleEqBl: return "Double Insert (equal_blocks)";
      }
    }

    /**
     * Try to insert one or two (cdag, c) operator pair(s) into the configuration
     *
     * @param config The configuration to update
     * @return On success, return the ratio between the inverse move probability
     *         and the move probability.
     *         On failure, return 0.0
     */
    double try_config_update(config_t &config) override;

    double insert_single_op_pair(long bl, config_t &config);

    scalar_t accept() override;

    fop_t last_d     = {};
    fop_t last_d_dag = {};

    acc_t &tau_diff_stat;
  };

} // namespace inchworm::moves

#include "insert_impl.hpp"
