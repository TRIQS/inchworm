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
  segment_t::segment_t(int pos1_, int pos2_, int id_) : pos1{pos1_}, pos2{pos2_}, id{id_}, size{pos2_ - pos1_} { EXPECTS(pos2 > pos1); }

  /**
   * Determine every possible segment within the diagram.
   * Each segment should:
   *  1. contain the same number of d and d_dag and
   *  2. not cross a split point
   *
   * Additionnal optimisation: a segment of length 4 and of type xoxo or oxox
   * does not need to be considered as it cannot be fully connected
   *
   * FIXME The full segment is currently always considered!
   */
  std::vector<segment_t> determine_segments(time_diagram_t const &diagram) {
    EXPECTS(not diagram.is_trivial); // FIXME Because of doublecounting of full segment?

    int diag_size = diagram.size();

    std::vector<segment_t> seg_list;

    int segment_id = 0;

    for (auto pos1 : range(0, diag_size - 1))                              //starting position of segment
      for (auto pos2 : range(pos1 + smallest_segment, diag_size + 1, 2)) { //one past the end of segment
	int len = pos2 - pos1;

	// No split-point should fall into the segment
        if (std::any_of(diagram.split_points.begin(), diagram.split_points.end(),
                        [&](auto &sp) { return pos1 < sp and sp < pos2; }))
          continue;

	// Optimization: Do not consider segments of length four
	if constexpr (remove_xoxo)
	  if (len == 4 and diagram.op_list[pos1].dag == diagram.op_list[pos1 + 2].dag) continue;

	// If [pos1, pos2) contains the same number of d and d_dag, add it to the list
	auto Ndag = std::count_if(diagram.op_list.cbegin() + pos1, diagram.op_list.cbegin() + pos2, [](auto &op) { return op.dag; });
	if (2 * Ndag == len) seg_list.emplace_back(pos1, pos2, segment_id++);
      }

    // FIXME WHY?
    //lastly, put the last segment (this one is k-connected and not fully connected. So we bypass the condition that it should not cross the split point:
    seg_list.emplace_back(0, diag_size, segment_id++);

    return seg_list;
  }

  set_of_segments_t::set_of_segments_t(segment_t const &seg0, time_diagram_t const &diagram)
     //FIXME set should be initialized empty
     : pos1{seg0.pos1}, pos2{seg0.pos2}, size{seg0.size}, segment_ids(diagram.perturbation_order() / (smallest_segment / 2)) {
    N_seg                = 0;
    segment_ids[N_seg++] = seg0.id;
  }

  void set_of_segments_t::append(segment_t const &seg1, time_diagram_t const &diagram) {
    //FIXME void set_of_segments_t::append(segment_t const &seg1) {
    EXPECTS(seg1.pos1 >= pos2);

    if (seg1.pos1 != pos2)
      adjacent = false;
    else if (std::any_of(begin(diagram.split_points), end(diagram.split_points), [j = pos2](int i) { return i == j; }))
      adjacent = false; // if the previous set of segments already ends at a split point, adding another one will make this set not adjacent anymore.
    else
      disjoint = false; // note, we consider that even if two segments touch at the split point, the set is still disjoint.

    size                 = seg1.pos2 - pos1;
    pos2                 = seg1.pos2;
    segment_ids[N_seg++] = seg1.id;
  }

  std::vector<set_of_segments_t> combine_segments(std::vector<segment_t> const &segment_list, time_diagram_t const &diagram, bool search_disjoint) {

    int n_seg = segment_list.size();

    std::vector<int> start_index_list = {0, 0};
    auto set_list                     = std::vector<set_of_segments_t>{}; // return value (pair[1])

    for (int j = 0; j < n_seg; j++) { set_list.push_back(set_of_segments_t(segment_list[j], diagram)); }
    for (int number_of_segment = 2; number_of_segment <= diagram.perturbation_order(); number_of_segment++) {
      start_index_list.push_back(set_list.size());

      int i1 = start_index_list.size() - 2;
      int i2 = start_index_list.size() - 1;
      for (int prev_index = start_index_list[i1]; prev_index < start_index_list[i2]; prev_index++) {

        set_of_segments_t previous_set = set_list[prev_index];
        //if we do not search for disjoint set, we search for adjacent set, only. We do not need the ones that are neither.

        for (auto const &additional_segment : segment_list) {
          // Only consider segments that start to the right of the current rightmost segment
          if (additional_segment.pos1 >= previous_set.pos2) {
            set_of_segments_t new_set = previous_set;
            new_set.append(additional_segment, diagram);

            if (search_disjoint) {
              if (not new_set.disjoint) continue;
            } else { // IMPORTANT distinction. A segment does not have to be disjoint or adjacent. But here, if we do not search for disjoint, we necessarly search for adjacent.
              if (not new_set.adjacent) continue;
            } //important brackets

            set_list.push_back(new_set);
          }
        }
      }
    }
    if (not search_disjoint) // if search adjacent erase the single segments (not necessary anymore)
      set_list.erase(set_list.begin(), set_list.begin() + start_index_list[2]);

    return set_list;
  }

  void calculate_segment(int segment_id,
                         std::vector<segment_t> &segment_list, // not const: modified
                         std::vector<set_of_segments_t> const &list_of_set_of_disjoint_segments,
                         std::vector<set_of_segments_t> const &list_of_set_of_adjacent_segments, hyb_matrix_t const &hyb_mat,
                         time_diagram_t const &diagram, bool special, int verbose) {

    auto &seg      = segment_list[segment_id];
    seg.calculated = true;
    if (verbose > 1) { print_segment(seg, diagram); }

    // [seg.pos1, ... , seg.pos2)
    sso_vector<int> range_of_vertex(seg.size);
    std::iota(range_of_vertex.begin(), range_of_vertex.end(), seg.pos1);

    if (verbose > 2) {
      std::printf("range of vertex: ");
      for (auto o : range_of_vertex) std::printf("%d ", o);
      hyb_mat.print();
    }

    // Calculate the determinant of the full segment
    seg.value += hyb_mat.extract_det(range_of_vertex);
    if (verbose > 2) std::printf("\nsegment[%d]= % 4.8f\n\n", segment_id, hyb_mat.extract_det(range_of_vertex));

    for (auto const &set : list_of_set_of_disjoint_segments) {
      if ((not special) and not((seg.pos1 <= set.pos1) and (seg.pos2 > set.pos2))) continue;

      if (special
          and (set.N_seg == 1) // this is the special case where we evaluate the full segment (at the end). We still need to exclude the itself.
          and ((seg.pos1 == set.pos1) and (seg.pos2 == set.pos2)))
        continue; // this is tricky, might have to change this at some point

      scalar_t value                   = 1.0;
      int sign_of_parcollet_charlebois = 1;
      bool is_finite                   = true;

      int number_of_vertex_to_remove = 0;
      for (int j = 0; j < set.N_seg; j++) number_of_vertex_to_remove += segment_list[set.segment_ids[j]].size;

      sso_vector<int> range_of_subvertex;
      range_of_subvertex.resize(seg.size - number_of_vertex_to_remove);

      int subseg_size = seg.size - number_of_vertex_to_remove;

      int index = 0;
      int pos   = range_of_vertex[0];

      for (int j = 0; j < set.N_seg; j++) {
        auto const &subseg = segment_list[set.segment_ids[j]];
        EXPECTS(seg.calculated);

        if (subseg.value == 0.0) { // somehow, this seems to happen often even if we consider float (does it still holds for complex numbers?)
          is_finite = false;
        }
        value *= -subseg.value;

        if (pos != subseg.pos1) {
          int number_of_new_vertex = (subseg.pos1 - pos);
          std::iota(range_of_subvertex.begin() + index, range_of_subvertex.begin() + index + number_of_new_vertex, pos);
          index += number_of_new_vertex;
        }
        pos = subseg.pos1 + subseg.size;

        if (subseg.size % 4 != 0)
          if ((subseg.pos2 - seg.pos1) % 2 == 1) sign_of_parcollet_charlebois *= -1;
      }
      std::iota(range_of_subvertex.begin() + index, range_of_subvertex.begin() + subseg_size, pos);

      if constexpr (remove_not_finite) {
        if (not is_finite) continue;
      } // avoid determinant calculation.

      if (range_of_subvertex.size() > 0) {
        scalar_t det1 = hyb_mat.extract_det(range_of_subvertex);
        value *= sign_of_parcollet_charlebois * det1;
      }
      seg.value += value;
    }

    seg.value_without_cuts = seg.value;

    for (auto const &cuts : list_of_set_of_adjacent_segments) {
      // Make sure that adjecent subset covers exactly the segment
      if (seg.pos1 == cuts.pos1)
        if (seg.pos2 == cuts.pos2) {

          scalar_t value = 1.0;
          for (int j = 0; j < cuts.N_seg; j++) {
            auto const &subseg = segment_list[cuts.segment_ids[j]];
            ASSERT(subseg.calculated);
            value *= -subseg.value_without_cuts;
          }
          seg.value -= value;
        }
    }
    if (verbose > 1) std::printf("   % 15.8f      % 15.8f\n", seg.value, seg.value_without_cuts);
  }

  scalar_t inclusion_exclusion(time_diagram_t const &diagram, hyb_matrix_t hyb_mat, int verbose) {

    if (diagram.size() == 0) return 1.0;

    //hyb_matrix_t hyb_mat(diagram, hyb_tau);
    hyb_mat.optimize_inclusion_exclusion(); // put some values to zero in hyb matrix (segment of length 2)
    if (diagram.is_trivial) { return 0; }   // not a question anymore: return 0 or det??
    if (diagram.perturbation_order() == 1) { return hyb_mat.det(); };

    if (verbose) std::printf("\n\n##################\nINCLUSION-EXCLUSION:\n");
    if (verbose > 1) {
      if (verbose > 2) hyb_mat.print();
    }
    std::vector<segment_t> segment_list = determine_segments(diagram);
    if (verbose) std::printf("\nsegment number = %lu\n\n", segment_list.size());

    std::vector<set_of_segments_t> list_of_set_of_disjoint_segments = combine_segments(segment_list, diagram, true);
    std::vector<set_of_segments_t> list_of_set_of_adjacent_segments = combine_segments(segment_list, diagram, false);

    if (verbose > 1) {
      std::printf("list of single segements:\n");
      print_diag(diagram);
      for (int j = 0; j < segment_list.size(); j++) {
        print_segment(segment_list[j], diagram);
        std::printf("\n");
      }
      std::printf("\n\nlist of set of disjoint segments:\n");
      print_diag(diagram);
      for (auto comb : list_of_set_of_disjoint_segments) {
        print_set(segment_list, comb, diagram);
        std::printf("\n");
      }
      std::printf("\n\nlist of set of adjacent segments:\n");
      print_diag(diagram);
      for (auto comb : list_of_set_of_adjacent_segments) {
        print_set(segment_list, comb, diagram);
        std::printf("\n");
      }
    }

    for (int length = smallest_segment; length <= 2 * diagram.perturbation_order(); length += 2) {
      for (auto seg : segment_list) {
        if (seg.size == length) {
          bool special = false;
          if (length == 2 * diagram.perturbation_order()) special = true;
          calculate_segment(seg.id, segment_list, list_of_set_of_disjoint_segments, list_of_set_of_adjacent_segments, hyb_mat, diagram, special,
                            verbose);
        }
      }
    }

    if (verbose > 0) {
      std::printf("\n## diagram = '%s'\n", diagram_string(diagram).c_str());
      std::printf("kOrder = %d\n", diagram.perturbation_order());
      std::printf("number of segments = %lu\n", segment_list.size());
      std::printf("number of adjacent sets = %lu\n", list_of_set_of_adjacent_segments.size());
      std::printf("number of disjoint sets = %lu\n", list_of_set_of_disjoint_segments.size());
    }

    return segment_list.back().value;
  }

  void print_segment(segment_t const &segment, time_diagram_t const &diagram) {
    std::vector<int> num_vector(diagram.size(), 0);
    for (int k = segment.pos1; k < segment.pos2; k++) num_vector[k] = 1;
    print_line(num_vector);
  }

  void print_set(std::vector<segment_t> const &segment_list, set_of_segments_t const &set_of_segments, time_diagram_t const &diagram) {
    std::vector<int> num_vector(diagram.size(), 0);

    for (int j = 0; j < set_of_segments.N_seg; j++) {
      for (int k = segment_list[set_of_segments.segment_ids[j]].pos1; k < segment_list[set_of_segments.segment_ids[j]].pos2; k++)
        num_vector[k] = j + 1;
    }
    print_line(num_vector);
  }

} // namespace inchworm::diagram
