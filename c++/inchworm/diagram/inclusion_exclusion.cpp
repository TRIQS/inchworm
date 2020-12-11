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
#include "./inclusion_exclusion.hpp"
#include "./utilities.hpp"
#include "./print.hpp"

namespace inchworm::diagram {
  segment_t::segment_t(int begin_, int end_, int id_) : begin{begin_}, end{end_}, id{id_}, size{end_ - begin_} { EXPECTS(end > begin); }

  set_of_segments_t::set_of_segments_t(segment_t const &seg, time_diagram_t const &diagram, std::vector<segment_t> const &segment_list)
     : begin{seg.begin},
       end{seg.end},
       size{seg.size},
       seg_ids_arr(diagram.perturbation_order() / (smallest_segment / 2)),
       N_segs{1},
       split_points_ptr{&diagram.split_points},
       segment_list_ptr{&segment_list} {
    seg_ids_arr[0] = seg.id;
  }

  void set_of_segments_t::append_right(segment_t const &seg) {
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

  std::vector<set_of_segments_t> combine_segments(std::vector<segment_t> const &segment_list, time_diagram_t const &diagram, bool search_disjoint) {

    auto set_list = std::vector<set_of_segments_t>{};
    for (auto const &seg : segment_list) { set_list.emplace_back(seg, diagram, segment_list); }
    if (not set_list.empty()) set_list.pop_back(); // Do not consider the full segment

    // Loop over the (continuously growing) list of sets until we have reached the end
    for (int set_idx = 0; set_idx != set_list.size(); ++set_idx) {

      for (auto const &additional_segment : segment_list) {

        // Only consider segments that start to the right of the current rightmost segment
        if (additional_segment.begin < set_list[set_idx].end) continue;

        set_of_segments_t new_set = set_list[set_idx];
        new_set.append_right(additional_segment);

        // Add set iff disjoint/adjacent when searching for disjoint/adjacent
        if ((search_disjoint && new_set.disjoint) or (not search_disjoint && new_set.adjacent)) { set_list.emplace_back(std::move(new_set)); }
      }
    }

    // If searching adjacent sets, erase all sets with a single segment
    if (not search_disjoint) set_list.erase(set_list.begin(), set_list.begin() + segment_list.size() - 1);

    return set_list;
  }

  void calculate_segment(int segment_id,
                         std::vector<segment_t> &segment_list, // not const: modified
                         std::vector<set_of_segments_t> const &list_of_set_of_disjoint_segments,
                         std::vector<set_of_segments_t> const &list_of_set_of_adjacent_segments, hyb_matrix_t const &hyb_mat,
                         time_diagram_t const &diagram, bool verbose) {
    auto &seg = segment_list[segment_id];
    EXPECTS(not seg.calculated);
    seg.calculated = true;

    // Calculate the determinant of the full segment
    seg.value = hyb_mat.extract_det(range(seg.begin, seg.end));

    // Vertex list for the det calculation. Allocate once and reuse
    std::vector<int> det_vertices;
    det_vertices.reserve(seg.size);

    for (auto const &set : list_of_set_of_disjoint_segments) {
      // FIXME Why do we have to make a distinction here between full segment and smaller? <= ?
      if ((seg.size < diagram.size()) and (set.begin < seg.begin or seg.end <= set.end)) continue;

      int pos        = seg.begin;
      scalar_t value = 1.0;

      for (auto const &subseg : set.segs()) {
        ASSERT(subseg.calculated);

        value *= -subseg.value;

        // Register any non-segment vertices before segment for det calculation
        for (auto i : range(pos, subseg.begin)) det_vertices.push_back(i);
        pos = subseg.end;

        // Caution: If pulling the sub-segment to the front of the current segment corresponds to
        // an odd number of row and column permutations in the segment submatrix we get an additional minus sign
        // Note: We could move this into the extract_det function
        if (subseg.size % 4 != 0 && (subseg.begin - seg.begin) % 2 == 1) value *= -1;
      }
      for (auto i : range(pos, seg.end)) det_vertices.push_back(i);

      if (det_vertices.size() > 0 && value != 0.0) { value *= hyb_mat.extract_det(det_vertices); }
      seg.value += value;

      det_vertices.clear();
    }

    seg.value_without_cuts = seg.value; // This should be the value in Eq (14) !?

    for (auto const &set : list_of_set_of_adjacent_segments) {
      if (seg.begin == set.begin && seg.end == set.end) {
        scalar_t value = 1.0;
        for (auto const &subseg : set.segs()) {
          ASSERT(subseg.calculated);
          value *= -subseg.value_without_cuts;
        }
        seg.value -= value;
      }
    }

    // ------------ Debug Prints ------------
    if (verbose) {
      print_segment(seg, diagram);
      std::printf("range of vertex: ");
      for (auto o : range(seg.begin, seg.end)) std::printf("%ld ", o);
      hyb_mat.print();
      std::printf("\nsegment[%d]= % 4.8f\n\n", segment_id, hyb_mat.extract_det(range(seg.begin, seg.end)));
      std::printf("   % 15.8f      % 15.8f\n", seg.value, seg.value_without_cuts);
    }
    // --------------------------------------
  }

  /**
   * Determine all relevant segment within the diagram.
   * This includes the full segment and each subsegment that
   *  1. contains the same number of d and d_dag operators
   *  2. does not cross a split point
   *
   * Additionnal optimisation: a segment of length 4 and of type xoxo or oxox
   * does not need to be considered as it cannot be fully connected
   */
  std::vector<segment_t> determine_segments(time_diagram_t const &diagram) {
    int diag_size = diagram.size();

    // Create sub-segments from smallest to largest, to respect calculation order
    int segment_id = 0;
    std::vector<segment_t> seg_list;
    for (auto size : range(smallest_segment, diag_size - 1, 2))
      for (auto begin : range(0, diag_size - size + 1)) {
        auto end = begin + size;

        // Optimization: Do not consider vanishing segments of size four (xoxo & oxox)
        if constexpr (remove_xoxo)
          if (size == 4 and diagram.op_list[begin].dag == diagram.op_list[begin + 2].dag) continue;

        // Consider only subsegments that do not contain a split-point
        if (std::any_of(diagram.split_points.begin(), diagram.split_points.end(), [&](auto &sp) { return begin < sp and sp < end; })) continue;

        // If [begin, end) contains the same number of d and d_dag, add it to the list
        auto Ndag = std::count_if(diagram.op_list.cbegin() + begin, diagram.op_list.cbegin() + end, [](auto &op) { return op.dag; });
        if (2 * Ndag == size) seg_list.emplace_back(begin, end, segment_id++);
      }

    // Finally, add the full segment, which may contain a split-point
    seg_list.emplace_back(0, diag_size, segment_id++);

    return seg_list;
  }

  scalar_t inclusion_exclusion(time_diagram_t const &diagram, hyb_matrix_t hyb_mat, bool verbose) {

    if (diagram.size() == 0) return 1.0;

    // Ignore segments of lenght 2 by setting necessary matrix elements to zero
    hyb_mat.optimize_inclusion_exclusion();

    // We need to have at least one split-point between operators for a finite hybridization weight
    if (diagram.is_trivial) { return 0; }

    if (diagram.perturbation_order() == 1) { return hyb_mat.det(); };

    std::vector<segment_t> segment_list                             = determine_segments(diagram);
    std::vector<set_of_segments_t> list_of_set_of_disjoint_segments = combine_segments(segment_list, diagram, true);
    std::vector<set_of_segments_t> list_of_set_of_adjacent_segments = combine_segments(segment_list, diagram, false);

    for (auto const &seg : segment_list) { // segment_list is sorted w.r.t. segment size
      calculate_segment(seg.id, segment_list, list_of_set_of_disjoint_segments, list_of_set_of_adjacent_segments, hyb_mat, diagram, verbose);
    }

    // ------------ Debug Prints ------------
    if (verbose) {
      std::printf("\n\n##################\nINCLUSION-EXCLUSION:\n");
      hyb_mat.print();
      std::printf("\nnumber of segments = %lu\n\n", segment_list.size());
      std::printf("list of single segements:\n");
      print_diag(diagram);
      for (int j = 0; j < segment_list.size(); j++) {
        print_segment(segment_list[j], diagram);
        std::printf("\n");
      }
      std::printf("\n\nlist of set of disjoint segments:\n");
      print_diag(diagram);
      for (auto const &s : list_of_set_of_disjoint_segments) {
        print_set(segment_list, s, diagram);
        std::printf("\n");
      }
      std::printf("\n\nlist of set of adjacent segments:\n");
      print_diag(diagram);
      for (auto const &s : list_of_set_of_adjacent_segments) {
        print_set(segment_list, s, diagram);
        std::printf("\n");
      }
      std::printf("\n## diagram = '%s'\n", diagram_string(diagram).c_str());
      std::printf("kOrder = %d\n", diagram.perturbation_order());
      std::printf("number of segments = %lu\n", segment_list.size());
      std::printf("number of adjacent sets = %lu\n", list_of_set_of_adjacent_segments.size());
      std::printf("number of disjoint sets = %lu\n", list_of_set_of_disjoint_segments.size());
    }
    // --------------------------------------

    return segment_list.back().value;
  }

  void print_segment(segment_t const &segment, time_diagram_t const &diagram) {
    std::vector<int> num_vector(diagram.size(), 0);
    for (auto k : range(segment.begin, segment.end)) num_vector[k] = 1;
    print_line(num_vector);
  }

  void print_set(std::vector<segment_t> const &segment_list, set_of_segments_t const &set, time_diagram_t const &diagram) {
    std::vector<int> num_vector(diagram.size(), 0);
    for (auto seg_id : set.seg_ids()) {
      auto const &seg = segment_list[seg_id];
      for (auto k : range(seg.begin, seg.end)) num_vector[k] = 1;
    }
    print_line(num_vector);
  }

} // namespace inchworm::diagram
