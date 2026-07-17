// Copyright (c) 2019--present, The Simons Foundation
// This file is part of inchworm and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once

#include "diagram.hpp"
#include "hyb_matrix.hpp"

namespace inchworm::diagram {

  /**
   * inclusion_exclusion algo based on Boag et al. PRB (2018) (with few changes)
   * We first determine the independent segments. We then combine
   * them into two lists: one fully disjoint (except for split points)
   * and another fully adjacent.
   * FIXME Extend documentation
   */
  hyb_scalar_t inclusion_exclusion(time_diagram_t const &diagram, hyb_matrix_t hyb_mat, bool verbose = false);

} // namespace inchworm::diagram
