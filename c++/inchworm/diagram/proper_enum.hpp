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

#include "diagram.hpp"
#include "hyb_matrix.hpp"

namespace inchworm::diagram {
  // find parity of a permutation by evaluating
  // by evaluating the parity of all cycle (or orbits)
  //
  int find_parity(std::vector<int> const &permutation);

  // print graph including diagram and all the arches for
  // a specific permutation. Different char are used to represent different status of the arch.
  //
  void print_graph(const std::vector<int> &permutation, const std::vector<bool> &cross_split_point, std::vector<bool> &visited,
                   time_diagram_t const &diagram);

  // Depth-first search (DFS) algorithm to search for every connected arch.
  // This recursive function will call itself until there is no more
  // free arch to visit (stored in variable visited).
  //
  void grow_pile(int arch, std::vector<int> const &permutation, std::vector<bool> &visited, time_diagram_t const &diagram);

  // DFS algorithm to search for every connected arch.
  // Everything is set up such that we can call the recursive DFS algorithm
  // "grow_pile" function. We need to first make a variable
  // that keep track of the arches that cross the split_point, we name
  // this vector cross_split_point_pile.
  //
  bool test_diagram_connection(const std::vector<int> &permutation, time_diagram_t const &diagram);

  // Look at all the arches combinaitions (represented by a permutation vector)
  // for a given diagram definition. We can use the function
  // "test_diagram_connection" to define if a diagram is proper or not.
  //
  scalar_t proper_enum(time_diagram_t const &diagram, hyb_matrix_t const &hyb_mat, int verbose = 0);

  // Just calculate the determinant using full enmuration
  //
  scalar_t full_enum(time_diagram_t const &diagram, hyb_matrix_t const &hyb_mat, int verbose = 0);

} // namespace inchworm::diagram
