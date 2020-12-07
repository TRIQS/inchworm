/*******************************************************************************
 *
 * inchworm: A TRIQS based impurity solver
 *
 * Copyright (c) 2019 The Simons foundation
 *   authors: Maxime Charlebois
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

#include "hyb_matrix.hpp"

namespace inchworm::diagram {

  /**
   * This class represents a single continuous segment with vertices in [pos1, pos2)
   *
   */
  struct segment_t {

    int pos1                    = 0;  // Index of the leftmost vertex
    int pos2                    = 0;  // Index one past the rightmost vertex
    int size                    = 0;  // The number of vertices in the segment
    scalar_t value              = 0.;
    scalar_t value_without_cuts = 0.;

    bool calculated = false; // Has this segment been calculated?
    int id          = 0;     // The segment identifier

    segment_t(int pos1_, int pos2_, int id_);
  };

  /**
   * Determine every possible segment based on the time_diagram definition.
   * The simple rule is: "Any segment should: 1. contain the same number of d
   * and d_dag and 2. not cross a split point".
   *
   * Additionnal optimisation: a segment of length 4 and of type xoxo or oxox
   * does not need to be considered as it cannot be fully connected
   */
  std::vector<segment_t> determine_segments(time_diagram_t const &diagram);

  /**
   * This class represents a set of non-overlapping segments contained in [pos1, pos2)
   *
   */
  struct set_of_segments_t {

    int pos1; // Index of the leftmost vertex

    int pos2;             // Index one past the rightmost vertex
    int size;             // Total length of the segment: pos2 - pos1
    bool disjoint = true; // The set is disjoint if all segments are separated by at least one vertex
    bool adjacent = true; // The set is adjacent if all segments touch and span the full range [pos1, pos2)

    sso_vector<int> segment_ids; // List that stores the subsegment indices. Initialize with maximal possible size.
    int N_seg;                   // Number of of segments, i.e. values in segment_ids that have been initialized

    // Constructor:
    //FIXME set_of_segments_t(time_diagram_t const &diagram);
    set_of_segments_t(segment_t const &seg0, time_diagram_t const &diagram);

    /**
     * Function to add a segment to the present set of segments.
     *
     * The function assumes that seg1 lies to the right of all segments contained so far.
     */
    void append(segment_t const &seg1, time_diagram_t const &diagram);
  };

  /**
   * This function generates the list of distjoint / adjoint sets
   * given the list of all possible segments
   *
   * It proceeds in steps. Every step reuses the previous set of segment. For exemple
   * when we try to generate set of 3 segments, we reuse every set
   * of 2 segment and try to append segments the segments in segment_list.
   * For this reason, we keep the information of the indices where the "N segments"
   * set start in the list "set_list". This is kept in the vector
   * start_index_list.
   */
  std::vector<set_of_segments_t> combine_segments(std::vector<segment_t> const &segment_list, time_diagram_t const &diagram, bool search_disjoint);

  /**
   * Calculate the value of one segment by analysing every segments of the set of segment (of both lists)
   *
   * The segment is coined special only if is the full segment of the diagram
   */
  void calculate_segment(int segment_id,
                         std::vector<segment_t> &segments_list, // not const: modified
                         std::vector<set_of_segments_t> const &list_of_set_of_disjoint_segments,
                         std::vector<set_of_segments_t> const &list_of_set_of_adjacent_segments, hyb_matrix_t const &hyb_mat,
                         time_diagram_t const &diagram, bool special = false, int verbose = 0);

  /**
   * inclusion_exclusion algo based on Boag et al. PRB (2018) (with few changes)
   * We first determine the independent segments. We then combine
   * them into two lists: one fully disjoint (except for split points)
   * and another fully adjacent.
   */
  scalar_t inclusion_exclusion(time_diagram_t const &diagram, hyb_matrix_t hyb_mat, int verbose = 0);

  /**
   * Print a single segment
   */
  void print_segment(segment_t const &segment, time_diagram_t const &diagram);

  /**
   * Print a set of segments
   */
  void print_set(std::vector<segment_t> const &segment_list, set_of_segments_t const &set_of_segments, time_diagram_t const &diagram);

} // namespace inchworm::diagram
