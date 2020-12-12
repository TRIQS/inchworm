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

#include <span>
#include <cppcoro/generator.hpp>

#include "./diagram.hpp"

namespace inchworm::diagram {

  // Global variables to set optimization levels
  constexpr int smallest_segment = 4;    // must be 2 or 4, beware.
  constexpr bool remove_xoxo     = true; // new optimisation 1: IMPORTANT, only works with the option smallest_segment = 4;

  /**
   * This class represents a single continuous segment with vertices in [begin, end)
   */
  struct segment_t {

    const int begin; // Index of the leftmost vertex
    const int end;   // Index one past the rightmost vertex
    const int id;    // The segment identifier
    const int size;  // The number of vertices in the segment

    bool calculated             = false; // Has this segment been calculated?
    scalar_t value              = 0.;    // The full value of the segment, c.f. Eq (14) arXiv:1807.00290v1
    scalar_t value_without_cuts = 0.;    // The value of the segment if split points are not considered // FIXME CHECK

    segment_t(int begin_, int end_, int id_) : begin{begin_}, end{end_}, id{id_}, size{end_ - begin_} { EXPECTS(end > begin); }
  };

  // -----------------------------------------------------------------------

  /**
   * This class represents a set of non-overlapping segments contained in [begin, end)
   */
  struct set_of_segments_t {

    int begin; // Index of the leftmost vertex
    int end;   // Index one past the rightmost vertex
    int size;  // Total length of the segment: end - begin

    bool disjoint = true; // The set is disjoint if all segments are separated by at least one vertex
    bool adjacent = true; // The set is adjacent if all segments touch and span the full diagram

    set_of_segments_t(segment_t const &seg, time_diagram_t const &diagram, std::vector<segment_t> const &segment_list)
       : begin{seg.begin},
         end{seg.end},
         size{seg.size},
         seg_ids_arr(diagram.perturbation_order() / (smallest_segment / 2)),
         N_segs{1},
         split_points_ptr{&diagram.split_points},
         segment_list_ptr{&segment_list} {
      seg_ids_arr[0] = seg.id;
    }

    std::span<int> seg_ids() { return {seg_ids_arr.data(), N_segs}; }
    std::span<const int> seg_ids() const { return {seg_ids_arr.data(), N_segs}; }

    cppcoro::generator<segment_t const &> segs() const {
      for (auto i : range(N_segs)) co_yield(*segment_list_ptr)[seg_ids_arr(i)];
    }

    /**
     * Function to add a segment to the present set of segments.
     * The function assumes that seg lies to the right of all segments contained so far.
     */
    void append_right(segment_t const &seg) {
      EXPECTS(seg.begin >= end);

      // We loose 'adjacency' if the new segment does not touch the previous right-most one or is separated by a split-point
      if (seg.begin > end or std::any_of(split_points_ptr->begin(), split_points_ptr->end(), [&](auto sp) { return sp == end; }))
        adjacent = false;
      else // We loose 'disjointness' if the new segment touches the previous right-most one and no split-point separates them
        disjoint = false;

      size                  = seg.end - begin;
      end                   = seg.end;
      seg_ids_arr[N_segs++] = seg.id;
    }

    private:
    sso_vector<int> seg_ids_arr; // Data array to store the subsegment indices. Initialized with maximal possible size.
    size_t N_segs;               // Number of of segments, i.e. values in data that have been initialized

    std::vector<int> const *split_points_ptr;       // Pointer to split_points vector of associated diagram
    std::vector<segment_t> const *segment_list_ptr; // Pointer to segment_list vector of associated diagram
  };

} // namespace inchworm::diagram
