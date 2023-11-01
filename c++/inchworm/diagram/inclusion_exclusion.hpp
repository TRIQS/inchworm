/*******************************************************************************
 *
 * inchworm: A TRIQS based impurity solver
 *
 * Copyright (c) 2019 The Simons foundation
 *   authors: Maxime Charlebois, Nils Wentzell
 *
 * inchworm is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * inchworm is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * inchworm. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
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
  // hyb_scalar_t inclusion_exclusion(time_diagram_t const &diagram, hyb_matrix_t hyb_mat, std::vector<segment_t>& segment_list,  std::vector<set_of_segments_t> &list_of_set_of_disjoint_segments, std::vector<set_of_segments_t> &list_of_set_of_adjacent_segments, bool verbose = false);

} // namespace inchworm::diagram
