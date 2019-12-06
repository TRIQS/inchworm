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

#include "./hybridization_function_matrix.hpp"
#include "./print.hpp"



// Definition of a segment:
//
struct segment_t {

  int pos1                                  = 0;
  int pos2                                  = 0; // note: by definition here, segment goes from index pos1 to pos2-1
  int size                                  = 0;
  hybridization_scalar_t value              = 0.;
  hybridization_scalar_t value_without_cuts = 0.;

  bool calculated = false;
  int numero      = 0;

  // Constructor:
  segment_t(int p1, int p2, int n) : pos1{p1}, pos2{p2}, numero{n} {
    EXPECTS(pos2 > pos1);
    size = pos2 - pos1;
  }
};


// print one segment 
// 
void print_segment(segment_t const &segment, time_diagram_t const &diagram) {
  std::vector<int> num_vector(diagram.list.size(), 0);
  for (int k = segment.pos1; k < segment.pos2; k++) num_vector[k] = 1;
  printLine(num_vector);
}


// Determine every possible segment based on the time_diagram definition.
// The simple rule is: "Any segment should: 1. contain the same number of c
// and cdag and 2. not cross a split point".
//
std::vector<segment_t> determine_segments(time_diagram_t const &diagram) {

  std::vector<segment_t> list;
  int N = diagram.list.size();

  int N_segment = 0;
  for (int i = 0; i < N - 1; i++)                           //starting position of segment
    for (int a = smallest_segment; a < N - i + 1; a += 2) { //length of segment

      if (std::any_of(diagram.split_points.begin(), diagram.split_points.end(),
                      [i, a](auto &split_point) { return segment_cross_p(i, i + a, split_point); }))
        continue;

      int Ndag = 0;
      for (int j = i; j < i + a; j++)
        if (diagram.list[j].dag) Ndag++;
      if (2 * Ndag == a) // check if same number of cdag an c in the segment starting at i and ending before i+a
      {
        auto seg = segment_t{i, i + a, N_segment++};
        list.push_back(seg);
        if (verbose > 0) {
          print_segment(seg, diagram);
          std::printf("\n");
        }
      }
    }

  //lastly, put the last segment (this one is k-connected and not fully connected. So we bypass the condition that it should not cross the split point:
  segment_t seg(0, N, N_segment++);
  list.push_back(seg);
  if(verbose) print_segment(seg, diagram);
  return list;
}


// Definition of combination_of_segments_t 
//
struct combination_of_segments_t {

  int pos1;
  int pos2; // note: by definition here, the combination of segments goes from index pos1 to pos2-1
  int size;
  bool disjoint = true; // we define disjoint when 2 segments does not touch (by convention, we choose a segment alone to be disjoint too)
  bool adjacent = true; // we define adjacent when all segments touches. If one does not, it is false.
  std::vector<int> list;
  
  // Constructor:
  combination_of_segments_t(segment_t const &seg0, time_diagram_t const &diagram) : pos1{seg0.pos1}, pos2{seg0.pos2}, size{seg0.size} {
    list.reserve(
       diagram
          .k_order() / (smallest_segment/2) ); //If smallest segment is length 2, we know that this is the maximum number of segments in a combination. If the smallest is 4, then it becomes k_order/2.
    list.push_back(seg0.numero);
  };

  // Function to add a segment to the present combination of segments:
  void append(segment_t const &seg1, time_diagram_t const &diagram) {
    if (seg1.pos1 != pos2)
      adjacent = false;
    else if (std::any_of(begin(diagram.split_points), end(diagram.split_points), [j = pos2](int i) { return i == j; }))
      adjacent =
         false; // if the previous combination of segments already ends at a split point, adding another one will make this combination not adjacent anymore.
    else
      disjoint = false; // note, we consider that even if two segments touch at the split point, the combination is still disjoint.

    size = seg1.pos2 - pos1;
    pos2 = seg1.pos2;
    list.push_back(seg1.numero);
  };
};


// print one combination of segments 
// 
void print_combin(std::vector<segment_t> const &segment_list, combination_of_segments_t const &combination_of_segments, time_diagram_t const &diagram) {
  std::vector<int> num_vector(diagram.list.size(), 0);

  for (int j = 0; j < combination_of_segments.list.size(); j++) {
    for (int k = segment_list[combination_of_segments.list[j]].pos1; 
              k < segment_list[combination_of_segments.list[j]].pos2; k++) num_vector[k] = j + 1;
  }
  printLine(num_vector);
}


// Combine the different segments defined in segment_list. It proceed in
// steps. Every step reuse the previous combination of segment. For exemple
// when we try to generate combination of 3 segments, we reuse every combination
// of 2 segment and try to append segments the segments in segment_list. 
// For this reason, we keep the information of the indices where the "N segments"
// combination start in the list "combination_list". This is kept in the vector
// start_index_list.
//
std::vector<combination_of_segments_t> combine_segments(std::vector<segment_t> const &segment_list, time_diagram_t const &diagram,
                                                        bool search_disjoint) {

  int n_seg = segment_list.size();

  std::vector<int> start_index_list = {0, 0};
  //std::vector<combination_of_segments_t> combination_list; // return value (pair[1])
  auto  combination_list = std::vector<combination_of_segments_t>{}; // return value (pair[1])
  

  for (int j = 0; j < n_seg; j++) { combination_list.push_back(combination_of_segments_t(segment_list[j], diagram)); }
  for (int number_of_segment = 2; number_of_segment <= diagram.k_order(); number_of_segment++) {
    start_index_list.push_back(combination_list.size());

    int i1 = start_index_list.size() - 2;
    int i2 = start_index_list.size() - 1;
    for (int prev_index = start_index_list[i1]; prev_index < start_index_list[i2]; prev_index++) {

      combination_of_segments_t previous_combination = combination_list[prev_index];
      //if we do not search for disjoint combination, we search for adjacent combination, only. We do not need the ones that are neither.

      for (auto additional_segment : segment_list) {
        if (additional_segment.pos1 >= previous_combination.pos2) {
          combination_of_segments_t new_combination = previous_combination;
          new_combination.append(additional_segment, diagram);

          if (search_disjoint) {
            if (not new_combination.disjoint) continue;
          } else { // IMPORTANT distinction. A segment does not have to be disjoint or adjacent. But here, if we do not search for disjoint, we necessarly search for adjacent. 
            if (not new_combination.adjacent) continue;
          } //important brackets

          combination_list.push_back(new_combination);
        }
      }
    }
  }
  return combination_list;
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
// by analysing every segiments of the combination of segment (of both lists)
//
void calculate_segment(int segment_numero,
                       std::vector<segment_t> &segments_list, // not const: modified
                       std::vector<combination_of_segments_t> const &combination_disjoint_list,
                       std::vector<combination_of_segments_t> const &combination_adjacent_list, hybridization_matrix const &hyb_mat, 
                       time_diagram_t const & diagram, 
                       bool special = false) {

  segments_list[segment_numero].calculated = true;
  if (verbose > 1) { print_segment(segments_list[segment_numero], diagram); }

  std::vector<int> range_of_vertex(segments_list[segment_numero].size);
  std::iota(range_of_vertex.begin(), range_of_vertex.end(), segments_list[segment_numero].pos1);

  segments_list[segment_numero].value += hyb_mat.extract_det(range_of_vertex);

  for (auto subs : combination_disjoint_list) {
    if ((not special) and not((segments_list[segment_numero].pos1 <= subs.pos1) and (segments_list[segment_numero].pos2 > subs.pos2))) continue;
    if (special and (subs.list.size() == 1)
        and ((segments_list[segment_numero].pos1 == subs.pos1) and (segments_list[segment_numero].pos2 == subs.pos2)))
      continue; // this is tricky, might have to change this at some point

    hybridization_scalar_t value = 1.0;
    std::vector<int> range_of_subvertex(range_of_vertex);
    int sign_of_parcollet_charlebois = 1;

    for (auto sub_segment_numero : subs.list) {
      segment_t seg = segments_list[sub_segment_numero];
      EXPECTS(seg.calculated);

      value *= -seg.value;

      range_of_subvertex = remove_segment_from_list(range_of_subvertex, seg.pos1, seg.size);

      if (seg.size % 4 != 0)
        if ((seg.pos2 - segments_list[segment_numero].pos1) % 2 == 1) sign_of_parcollet_charlebois *= -1;
    }

    if (range_of_subvertex.size() > 0) {
      hybridization_scalar_t det1 = hyb_mat.extract_det(range_of_subvertex);
      value *= sign_of_parcollet_charlebois * det1;
    }
    segments_list[segment_numero].value += value;
  }

  segments_list[segment_numero].value_without_cuts = segments_list[segment_numero].value;

  for (auto cuts : combination_adjacent_list) {
    if (segments_list[segment_numero].pos1 == cuts.pos1)
      if (segments_list[segment_numero].pos2 == cuts.pos2)
        if (cuts.list.size() > 1) {

          hybridization_scalar_t value = 1.0;

          for (auto sub_segment_numero : cuts.list) {
            segment_t seg = segments_list[sub_segment_numero];
            EXPECTS(seg.calculated);

            value *= -seg.value_without_cuts;
          }
          segments_list[segment_numero].value -= value;
        }
  }
  if (verbose > 1) std::printf("   % 15.8f      % 15.8f\n", segments_list[segment_numero].value, segments_list[segment_numero].value_without_cuts);
}


// inclusion_exclusion algo based on Boag et al. PRB (2018) (with few changes)
// We first determine the independent segments. We then combine
// them into two lists: one fully disjoint (except for split points) 
// and another fully adjacent.
//
hybridization_scalar_t inclusion_exclusion(time_diagram_t const &diagram) {

  if (verbose) print_diag(diagram);
  std::vector<segment_t> segment_list = determine_segments(diagram);
  if (verbose) std::printf("\nsegment number = %lu\n\n", segment_list.size());

  std::vector<combination_of_segments_t> combination_disjoint_list = combine_segments(segment_list, diagram, true);
  std::vector<combination_of_segments_t> combination_adjacent_list = combine_segments(segment_list, diagram, false);

  if (verbose > 1) {
    std::printf("\n\ncombination of disjoint segments:\n\n");
    print_diag(diagram);
    for (auto comb : combination_disjoint_list) {
      print_combin(segment_list, comb, diagram);
      std::printf("\n");
    }
    std::printf("\n\ncombination of adjacent segments:\n\n");
    print_diag(diagram);
    for (auto comb : combination_adjacent_list) {
      print_combin(segment_list, comb, diagram);
      std::printf("\n");
    }
  }

  hybridization_matrix hyb_mat(diagram);

  for (int length = smallest_segment; length <= 2 * diagram.k_order(); length += 2) {
    if (verbose > 1) std::printf("\n############\nsegment length = %d\n", length);
    for (auto seg : segment_list) {
      if (seg.size == length) {
        bool special = false;
        if (length == 2 * diagram.k_order()) special = true;
        calculate_segment(seg.numero, segment_list, combination_disjoint_list, combination_adjacent_list, hyb_mat, diagram, special);
      }
    }
  }

  if (verbose > 0) {
    std::printf("\n## diagram = '%s'\n", diagram_string(diagram).c_str());
    std::printf("kOrder = %d\n", diagram.k_order());
    std::printf("number of segments = %lu\n", segment_list.size());
    std::printf("number of adjacent combinations = %lu\n", combination_adjacent_list.size());
    std::printf("number of disjoint combinations = %lu\n", combination_disjoint_list.size());
  }

  return segment_list.back().value;
}

