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
#include "./hyb_matrix.hpp"
#include "./print.hpp"

namespace inchworm::diagram {
  // Definition of a segment:
  //
  struct segment_t {

    int pos1                    = 0;
    int pos2                    = 0; // note: by definition here, segment goes from index pos1 to pos2-1
    int size                    = 0;
    scalar_t value              = 0.;
    scalar_t value_without_cuts = 0.;

    bool calculated = false;
    int numero      = 0;

    // Constructor:
    segment_t(int p1, int p2, int n);
  };

  // print one segment
  //
  void print_segment(segment_t const &segment, time_diagram_t const &diagram);

  // Determine every possible segment based on the time_diagram definition.
  // The simple rule is: "Any segment should: 1. contain the same number of c
  // and d_dag and 2. not cross a split point".
  //
  // Additionnal optimisation: a segment of length 4 and of type xoxo or oxox
  // does not need to be considered as it cannot be fully connected
  //
  std::vector<segment_t> determine_segments(time_diagram_t const &diagram);

  // Definition of set_of_segments_t
  //
  struct set_of_segments_t {

    int pos1;
    int pos2; // note: by definition here, the set of segments goes from index pos1 to pos2-1
    int size;
    bool disjoint = true; // we define disjoint when 2 segments does not touch (by convention, we choose a segment alone to be disjoint too)
    bool adjacent = true; // we define adjacent when all segments touches. If one does not, it is false.
    std::vector<int> list;

    // Constructor:
    set_of_segments_t(segment_t const &seg0, time_diagram_t const &diagram);

    // Function to add a segment to the present set of segments:
    void append(segment_t const &seg1, time_diagram_t const &diagram);
  };

  // print one set of segments
  //
  void print_set(std::vector<segment_t> const &segment_list, set_of_segments_t const &set_of_segments, time_diagram_t const &diagram);

  // Combine the different segments defined in segment_list. It proceed in
  // steps. Every step reuse the previous set of segment. For exemple
  // when we try to generate set of 3 segments, we reuse every set
  // of 2 segment and try to append segments the segments in segment_list.
  // For this reason, we keep the information of the indices where the "N segments"
  // set start in the list "set_list". This is kept in the vector
  // start_index_list.
  //
  std::vector<set_of_segments_t> combine_segments(std::vector<segment_t> const &segment_list, time_diagram_t const &diagram, bool search_disjoint);

  // Remove "segment_size" elements of a vector of int, starting at the value "segment_min".
  // Of course is requires the "segment_min" to be in the vector and more than "segment_size"
  // away from the end of the vector. This is ensured by the EXPECTS() check.
  //
  std::vector<int> remove_segment_from_list(std::vector<int> const &list, int segment_min, int segment_size);

  // Calculate the value of one segment
  // by analysing every segiments of the set of segment (of both lists)
  //
  void calculate_segment(int segment_numero,
                         std::vector<segment_t> &segments_list, // not const: modified
                         std::vector<set_of_segments_t> const &set_disjoint_list, std::vector<set_of_segments_t> const &set_adjacent_list,
                         hyb_matrix_t const &hyb_mat, time_diagram_t const &diagram, bool special = false, int verbose = 0);

  // inclusion_exclusion algo based on Boag et al. PRB (2018) (with few changes)
  // We first determine the independent segments. We then combine
  // them into two lists: one fully disjoint (except for split points)
  // and another fully adjacent.
  //
  //scalar_t inclusion_exclusion(time_diagram_t const &diagram, hyb_adaptor_t const &hyb_tau);
  //scalar_t inclusion_exclusion(time_diagram_t const &diagram, std::function<scalar_t(double)>);

  //scalar_t determinant(time_diagram_t const &diagram, hyb_adaptor_t const &hyb_tau);
  //scalar_t determinant(time_diagram_t const &diagram, std::function<scalar_t(double)>);
  scalar_t inclusion_exclusion(time_diagram_t const &diagram, hyb_matrix_t hyb_mat, int verbose = 0);
} // namespace inchworm::diagram
