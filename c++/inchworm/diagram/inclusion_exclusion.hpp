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
   * This class represents a single continuous segment with vertices in [begin, end)
   *
   */
  struct segment_t {

    const int begin; // Index of the leftmost vertex
    const int end;   // Index one past the rightmost vertex
    const int id;    // The segment identifier
    const int size;  // The number of vertices in the segment

    bool calculated             = false; // Has this segment been calculated?
    scalar_t value              = 0.;    // The full value of the segment, c.f. Eq (14) arXiv:1807.00290v1
    scalar_t value_without_cuts = 0.;    // The value of the segment if split points are not considered // FIXME CHECK

    segment_t(int begin_, int end_, int id_);
  };

  /**
   * This class represents a set of non-overlapping segments contained in [begin, end)
   *
   */
  struct set_of_segments_t {

    int begin; // Index of the leftmost vertex
    int end;   // Index one past the rightmost vertex
    int size;  // Total length of the segment: end - begin

    bool disjoint = true; // The set is disjoint if all segments are separated by at least one vertex
    bool adjacent = true; // The set is adjacent if all segments touch and span the full diagram

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
                         time_diagram_t const &diagram, bool special = false, bool verbose = false);

  /**
   * inclusion_exclusion algo based on Boag et al. PRB (2018) (with few changes)
   * We first determine the independent segments. We then combine
   * them into two lists: one fully disjoint (except for split points)
   * and another fully adjacent.
   */
  scalar_t inclusion_exclusion(time_diagram_t const &diagram, hyb_matrix_t hyb_mat, bool verbose = false);

  /**
   * Print a single segment
   */
  void print_segment(segment_t const &segment, time_diagram_t const &diagram);

  /**
   * Print a set of segments
   */
  void print_set(std::vector<segment_t> const &segment_list, set_of_segments_t const &set_of_segments, time_diagram_t const &diagram);

} // namespace inchworm::diagram
