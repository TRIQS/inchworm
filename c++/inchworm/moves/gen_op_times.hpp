// Copyright (c) 2021--present, The Simons Foundation
// This file is part of inchworm and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once

#include "./../config.hpp"

#include <triqs/mc_tools/random_generator.hpp>

namespace inchworm::moves {

  std::array<double, 2> gen_op_times(triqs::mc_tools::random_generator &rng, fop_t const &d, fop_t const &d_dag, config_t const &config,
                                     double tau_split, double tau_max);

  double get_time_prop_prob(fop_t const &d, fop_t const &d_dag, config_t const &config, double tau_split, double tau_max);

} // namespace inchworm::moves
