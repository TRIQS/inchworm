// Copyright (c) 2020--present, The Simons Foundation
// This file is part of inchworm and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once

#include "base_move.hpp"

namespace inchworm::moves {

  template <KIND Kind> struct remove : public base_move {
    using base_move::base_move;

    inline std::string name() const override {
      switch (Kind) {
        case KIND::Single: return "Remove";
        case KIND::Double: return "Double Remove";
        case KIND::DoubleEqBl: return "Double Remove (equal_blocks)";
      }
    }

    /**
     * Try to remove two or one (cdag, c) operator pair(s) from the configuration
     *
     * @param config The configuration to update
     * @return On success, return the ratio between the inverse move probability
     *         and the move probability.
     *         On failure, return 0.0
     */
    double try_config_update(config_t &config) override;

    double remove_single_op_pair(long bl, config_t &config);
  };

} // namespace inchworm::moves

#include "remove_impl.hpp"
