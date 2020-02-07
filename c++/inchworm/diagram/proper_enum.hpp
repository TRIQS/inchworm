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
#include <iostream>
#include <string>
#include <utility>

#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <assert.h> // will be removed later
#include <vector>

#include "diagram.hpp"
#include "utilities.hpp"
#include "print.hpp"
#include "hyb_matrix.hpp"

namespace inchworm::diagram {
  // find parity of a permutation by evaluating
  // by evaluating the parity of all cycle (or orbits)
  //
  int find_parity(std::vector<int> const &permutation) {
    int parity = 1;
    std::vector<bool> visited(permutation.size(), false);

    for (int ii = 0; ii < permutation.size(); ii++) {
      if (visited[ii]) continue;
      int jj = permutation[ii];
      while (jj != ii) {
        parity *= -1; //change for addiiton?
        visited[jj] = true;
        jj          = permutation[jj];
      }
    }
    return parity;
  }

  /*
// return a vector of position where the char "c" appear in the string "str"
//
std::vector<int> find_positions(std::string const &str, char c) {
  std::vector<int> v;
  for (int i = 0; i < str.size(); i++)
    if (str[i] == c) v.push_back(i);
  return v;
}
*/

  // print graph including diagram and all the arches for
  // a specific permutation. Different char are used to represent different status of the arch.
  //
  void print_graph(const std::vector<int> &permutation, const std::vector<bool> &cross_split_point, std::vector<bool> &visited,
                   time_diagram_t const &diagram) {

    print_diag(diagram);

    int arch;
    for (arch = 0; arch < diagram.perturbation_order(); arch++) { //as much arches than pair of vertices
      char char1 = '.';
      if (cross_split_point[arch]) char1 = '_';
      if (visited[arch]) char1 = '=';
      printArch(diagram.pos_c[arch], diagram.pos_cdag[permutation[arch]], diagram.perturbation_order(), char1);
    }
    std::printf("\n");
  }

  // Deep First Search (DFS) algorithm to search for every connected arch.
  // This recursive function will call itself until there is no more
  // free arch to visit (stored in variable visited).
  //
  void grow_pile(int arch, std::vector<int> const &permutation, std::vector<bool> &visited, time_diagram_t const &diagram) {

    //connexion_pile.push_back(arch);
    int c    = diagram.pos_c[arch];
    int cdag = diagram.pos_cdag[permutation[arch]];

    for (int new_arch = 0; new_arch < diagram.perturbation_order(); new_arch++) {
      if (visited[new_arch]) continue;

      if (segment_cross(c, cdag, diagram.pos_c[new_arch], diagram.pos_cdag[permutation[new_arch]])) {
        visited[new_arch] = true;
        grow_pile(new_arch, permutation, visited, diagram);
      }
    }
    //connexion_pile.pop_back();
    visited[arch] = true;
  }

  // DFS algorithm to search for every connected arch.
  // Everything is set up such that we can call the recursive DFS algorithm
  // "grow_pile" function. We need to first make a variable
  // that keep track of the arches that cross the split_point, we name
  // this vector cross_split_point_pile.
  //
  bool test_diagram_connection(const std::vector<int> &permutation, time_diagram_t const &diagram) {

    std::vector<bool> cross_split_point(diagram.perturbation_order(), false); // for graphic purpose only (to change)
    std::vector<bool> visited(diagram.perturbation_order(), false);

    std::vector<int> cross_split_point_pile;
    cross_split_point_pile.reserve(diagram.perturbation_order()); // we know that the pile will not grow bigger than the number of arch = k_order

    for (int arch = 0; arch < diagram.perturbation_order(); arch++) {
      int c    = diagram.pos_c[arch];
      int cdag = diagram.pos_cdag[permutation[arch]];

      //for(auto p : diagram.split_points) std::printf("c=%d  cdag=%d   p=%d  \n", c, cdag, p);

      if (std::any_of(diagram.split_points.begin(), diagram.split_points.end(),
                      [c, cdag](auto &split_point) { return segment_cross_point(c, cdag, split_point - 1); })) {
        //if (segment_cross_point(c, cdag, diagram.split_points)) {
        cross_split_point[arch] = true;
        cross_split_point_pile.push_back(arch);
      }
    }

    if constexpr (verbose > 2) {
      std::printf("\n");
      print_graph(permutation, cross_split_point, visited, diagram);
    }

    for (auto arch : cross_split_point_pile) {
      if (visited[arch]) continue;
      visited[arch] = true;
      grow_pile(arch, permutation, visited, diagram);
    }

    if constexpr (verbose > 2) print_graph(permutation, cross_split_point, visited, diagram);

    return std::all_of(visited.begin(), visited.end(), [](bool v) { return v; });
  }

  // Look at all the arches combinaitions (represented by a permutation vector)
  // for a given diagram definition. We can use the function
  // "test_diagram_connection" to define if a diagram is proper or not.
  //
  scalar_t proper_enum(time_diagram_t const &diagram, hyb_matrix_t const &hyb_mat) {

    if constexpr (verbose) std::printf("\n\n##################\nPROPER-ENUMERATION:\n");
    //auto hyb_mat     = hyb_matrix_t{diagram, hyb_function};
    auto permutation = std::vector<int>(diagram.perturbation_order());

    for (int i = 0; i < diagram.perturbation_order(); i++) permutation[i] = i;

    int NN = 0, N_proper = 0;
    scalar_t total_value = 0.0, value = 1.0;
    do {
      NN += 1;
      int parity = find_parity(permutation);
      value      = 1.0;
      for (int i = 0; i < permutation.size(); i++) value *= hyb_mat.mat(permutation[i], i);

      if constexpr (verbose > 2) {
        std::printf("\ndiagram #%d:  permutation (", NN);
        for (auto i : permutation) std::printf("%d ", i);
        std::printf(")");
      }
      if (test_diagram_connection(permutation, diagram)) {
        N_proper += 1;
        total_value += parity * value;
        if constexpr (verbose > 2) std::printf("proper,   value=% 4.7f, parity=%d\n\n", value, parity);
      } else if constexpr (verbose > 2)
        std::printf("improper, value=% 4.7f, parity=%d\n\n", value, parity);
    } while (std::next_permutation(permutation.begin(), permutation.end()));

    if constexpr (verbose) {
      std::printf("## diagram = '%s'\n", diagram_string(diagram).c_str());
      std::printf("k_order = %d\n", diagram.perturbation_order());
      std::printf("number of diagram = %d\n", NN);
      std::printf("number of proper diagram = %d\n", N_proper);
    }
    return total_value;
  }
} // namespace inchworm::diagram
