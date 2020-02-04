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

namespace inchworm::diagram {
  segment_t::segment_t(int p1, int p2, int n) : pos1{p1}, pos2{p2}, numero{n} {
    EXPECTS(pos2 > pos1);
    size = pos2 - pos1;
  }

  // print one segment
  //
  void print_segment(segment_t const &segment, time_diagram_t const &diagram) {
    std::vector<int> num_vector(diagram.op_list.size(), 0);
    for (int k = segment.pos1; k < segment.pos2; k++) num_vector[k] = 1;
    printLine(num_vector);
  }

  // Determine every possible segment based on the time_diagram definition.
  // The simple rule is: "Any segment should: 1. contain the same number of c
  // and cdag and 2. not cross a split point".
  //
  // Additionnal optimisation: a segment of length 4 and of type xoxo or oxox
  // does not need to be considered as it cannot be fully connected
  //
  std::vector<segment_t> determine_segments(time_diagram_t const &diagram) {

    std::vector<segment_t> seg_list;
    int N = diagram.op_list.size();

    int N_segment = 0;
    for (int i = 0; i < N - 1; i++)                           //starting position of segment
      for (int a = smallest_segment; a < N - i + 1; a += 2) { //length of segment

        if (std::any_of(diagram.split_points.begin(), diagram.split_points.end(),
                        [i, a](auto &split_point) { return segment_cross_p(i, i + a, split_point); }))
          continue;

        int Ndag = 0;

        if constexpr (remove_xoxo)
          if (a == 4 and diagram.op_list[i].dag == diagram.op_list[i + 2].dag) continue;

        // count the number of dag, must be half of the lenght a:
        for (int j = i; j < i + a; j++)
          if (diagram.op_list[j].dag) Ndag++;
        if (2 * Ndag == a) // check if same number of cdag an c in the segment starting at i and ending before i+a
        {
          auto seg = segment_t{i, i + a, N_segment++};
          seg_list.push_back(seg);
          if constexpr (verbose > 0) {
            print_segment(seg, diagram);
            std::printf("\n");
          }
        }
      }

    //lastly, put the last segment (this one is k-connected and not fully connected. So we bypass the condition that it should not cross the split point:
    segment_t seg(0, N, N_segment++);
    seg_list.push_back(seg);
    if constexpr (verbose) print_segment(seg, diagram);
    return seg_list;
  }

  set_of_segments_t::set_of_segments_t(segment_t const &seg0, time_diagram_t const &diagram) : pos1{seg0.pos1}, pos2{seg0.pos2}, size{seg0.size} {
    list.reserve(
       diagram.perturbation_order()
       / (smallest_segment
          / 2)); //If smallest segment is length 2, we know that this is the maximum number of segments in a set. If the smallest is 4, then it becomes k_order/2.
    list.push_back(seg0.numero);
  }

  // Function to add a segment to the present set of segments:
  void set_of_segments_t::append(segment_t const &seg1, time_diagram_t const &diagram) {
    if (seg1.pos1 != pos2)
      adjacent = false;
    else if (std::any_of(begin(diagram.split_points), end(diagram.split_points), [j = pos2](int i) { return i == j; }))
      adjacent = false; // if the previous set of segments already ends at a split point, adding another one will make this set not adjacent anymore.
    else
      disjoint = false; // note, we consider that even if two segments touch at the split point, the set is still disjoint.

    size = seg1.pos2 - pos1;
    pos2 = seg1.pos2;
    list.push_back(seg1.numero);
  }

  // print one set of segments
  //
  void print_set(std::vector<segment_t> const &segment_list, set_of_segments_t const &set_of_segments, time_diagram_t const &diagram) {
    std::vector<int> num_vector(diagram.op_list.size(), 0);

    for (int j = 0; j < set_of_segments.list.size(); j++) {
      for (int k = segment_list[set_of_segments.list[j]].pos1; k < segment_list[set_of_segments.list[j]].pos2; k++) num_vector[k] = j + 1;
    }
    printLine(num_vector);
  }

  // Combine the different segments defined in segment_list. It proceed in
  // steps. Every step reuse the previous set of segment. For exemple
  // when we try to generate set of 3 segments, we reuse every set
  // of 2 segment and try to append segments the segments in segment_list.
  // For this reason, we keep the information of the indices where the "N segments"
  // set start in the list "set_list". This is kept in the vector
  // start_index_list.
  //
  std::vector<set_of_segments_t> combine_segments(std::vector<segment_t> const &segment_list, time_diagram_t const &diagram, bool search_disjoint) {

    int n_seg = segment_list.size();

    std::vector<int> start_index_list = {0, 0};
    //std::vector<set_of_segments_t> set_list; // return value (pair[1])
    auto set_list = std::vector<set_of_segments_t>{}; // return value (pair[1])

    for (int j = 0; j < n_seg; j++) { set_list.push_back(set_of_segments_t(segment_list[j], diagram)); }
    for (int number_of_segment = 2; number_of_segment <= diagram.perturbation_order(); number_of_segment++) {
      start_index_list.push_back(set_list.size());

      int i1 = start_index_list.size() - 2;
      int i2 = start_index_list.size() - 1;
      for (int prev_index = start_index_list[i1]; prev_index < start_index_list[i2]; prev_index++) {

        set_of_segments_t previous_set = set_list[prev_index];
        //if we do not search for disjoint set, we search for adjacent set, only. We do not need the ones that are neither.

        for (auto additional_segment : segment_list) {
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

  // Remove "segment_size" elements of a vector of int, starting at the value "segment_min".
  // Of course is requires the "segment_min" to be in the vector and more than "segment_size"
  // away from the end of the vector. This is ensured by the EXPECTS() check.
  //
  std::vector<int> remove_segment_from_list(std::vector<int> const &list, int segment_min, int segment_size) {
    std::vector<int> output;
    output.reserve(list.size() - segment_size);

    int segment_max = segment_min + segment_size;
    for (auto l : list)
      if ((l < segment_min) or (l >= segment_max)) { output.push_back(l); }

    EXPECTS(output.size() == list.size() - segment_size)
    return output;
  }

  // Calculate the value of one segment
  // by analysing every segiments of the set of segment (of both lists)
  //
  void calculate_segment(int segment_numero,
                         std::vector<segment_t> &segments_list, // not const: modified
                         std::vector<set_of_segments_t> const &set_disjoint_list, std::vector<set_of_segments_t> const &set_adjacent_list,
                         hyb_matrix_t const &hyb_mat, time_diagram_t const &diagram, bool special) {

    segments_list[segment_numero].calculated = true;
    if constexpr (verbose > 1) { print_segment(segments_list[segment_numero], diagram); }

    std::vector<int> range_of_vertex(segments_list[segment_numero].size);
    std::iota(range_of_vertex.begin(), range_of_vertex.end(), segments_list[segment_numero].pos1);

    segments_list[segment_numero].value += hyb_mat.extract_det(range_of_vertex);

    for (auto subs : set_disjoint_list) {
      if ((not special) and not((segments_list[segment_numero].pos1 <= subs.pos1) and (segments_list[segment_numero].pos2 > subs.pos2))) continue;
      if (special and (subs.list.size() == 1)
          and ((segments_list[segment_numero].pos1 == subs.pos1) and (segments_list[segment_numero].pos2 == subs.pos2)))
        continue; // this is tricky, might have to change this at some point

      scalar_t value = 1.0;
      std::vector<int> range_of_subvertex(range_of_vertex);
      int sign_of_parcollet_charlebois = 1;
      bool is_finite                   = true;

      for (auto sub_segment_numero : subs.list) {
        segment_t seg = segments_list[sub_segment_numero];
        EXPECTS(seg.calculated);

        if (seg.value == 0.0) { // somehow, this seems to happen often even if we consider float (does it still holds for complex numbers?)
          is_finite = false;
          //std::printf("is not finite: %e \n", seg.value);
        }
        value *= -seg.value;

        range_of_subvertex = remove_segment_from_list(range_of_subvertex, seg.pos1, seg.size);

        if (seg.size % 4 != 0)
          if ((seg.pos2 - segments_list[segment_numero].pos1) % 2 == 1) sign_of_parcollet_charlebois *= -1;
      }

      if constexpr (remove_not_finite) {
        if (not is_finite) continue;
      } // avoid determinant calculation.

      if (range_of_subvertex.size() > 0) {
        scalar_t det1 = hyb_mat.extract_det(range_of_subvertex);
        value *= sign_of_parcollet_charlebois * det1;
      }
      segments_list[segment_numero].value += value;
    }

    segments_list[segment_numero].value_without_cuts = segments_list[segment_numero].value;

    for (auto cuts : set_adjacent_list) {
      if (segments_list[segment_numero].pos1 == cuts.pos1)
        if (segments_list[segment_numero].pos2 == cuts.pos2) {

          scalar_t value = 1.0;

          for (auto sub_segment_numero : cuts.list) {
            segment_t seg = segments_list[sub_segment_numero];
            EXPECTS(seg.calculated);

            value *= -seg.value_without_cuts;
          }
          segments_list[segment_numero].value -= value;
        }
    }
    if constexpr (verbose > 1)
      std::printf("   % 15.8f      % 15.8f\n", segments_list[segment_numero].value, segments_list[segment_numero].value_without_cuts);
  }

  // inclusion_exclusion algo based on Boag et al. PRB (2018) (with few changes)
  // We first determine the independent segments. We then combine
  // them into two lists: one fully disjoint (except for split points)
  // and another fully adjacent.
  //
  scalar_t inclusion_exclusion(time_diagram_t const &diagram, std::function<scalar_t(double)> hyb_function) {

    hyb_matrix_t hyb_mat(diagram, hyb_function);
    hyb_mat.optimize_inclusion_exclusion(); // put some values to zero in hyb matrix (segment of length 2)
    if (diagram.is_trivial) { return 0; }   // return 0 or det??
    if (diagram.perturbation_order() == 1) {
      if (std::any_of(begin(diagram.split_points), end(diagram.split_points), [](int i) { return i == 1; }))
        return hyb_mat.det();
      else
        return 0.0;
    };

    if constexpr (verbose) {
      std::printf("\n\n##################\nINCLUSION-EXCLUSION:\n");
      std::printf("list of single segements:\n");
      print_diag(diagram);
    }
    std::vector<segment_t> segment_list = determine_segments(diagram);
    if constexpr (verbose) std::printf("\nsegment number = %lu\n\n", segment_list.size());

    std::vector<set_of_segments_t> set_disjoint_list = combine_segments(segment_list, diagram, true);
    std::vector<set_of_segments_t> set_adjacent_list = combine_segments(segment_list, diagram, false);

    if constexpr (verbose > 1) {
      std::printf("\n\nlist of set of disjoint segments:\n");
      print_diag(diagram);
      for (auto comb : set_disjoint_list) {
        print_set(segment_list, comb, diagram);
        std::printf("\n");
      }
      std::printf("\n\nlist of set of adjacent segments:\n");
      print_diag(diagram);
      for (auto comb : set_adjacent_list) {
        print_set(segment_list, comb, diagram);
        std::printf("\n");
      }
    }

    for (int length = smallest_segment; length <= 2 * diagram.perturbation_order(); length += 2) {
      if constexpr (verbose > 1) std::printf("\n############\nsegment length = %d\n", length);
      for (auto seg : segment_list) {
        if (seg.size == length) {
          bool special = false;
          if (length == 2 * diagram.perturbation_order()) special = true;
          calculate_segment(seg.numero, segment_list, set_disjoint_list, set_adjacent_list, hyb_mat, diagram, special);
        }
      }
    }

    if constexpr (verbose > 0) {
      std::printf("\n## diagram = '%s'\n", diagram_string(diagram).c_str());
      std::printf("kOrder = %d\n", diagram.perturbation_order());
      std::printf("number of segments = %lu\n", segment_list.size());
      std::printf("number of adjacent sets = %lu\n", set_adjacent_list.size());
      std::printf("number of disjoint sets = %lu\n", set_disjoint_list.size());
    }

    return segment_list.back().value;
  }

  scalar_t determinant(time_diagram_t const &diagram, std::function<scalar_t(double)> hyb_function) {
    hyb_matrix_t hyb_mat(diagram, hyb_function);
    //hyb_mat.print();
    return hyb_mat.det();
  }

} // namespace inchworm::diagram
