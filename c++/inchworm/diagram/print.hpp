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

#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <numeric>

#include "./diagram.hpp"

///////////////////////////////////
////////// PRINT //////////////////
///////////////////////////////////

namespace inchworm::diagram {
  void print_vector(std::vector<int> const &v);

  // print diagram and its split point above.
  //
  void print_diag(time_diagram_t const &diagram);

  void print_configuration(time_diagram_t const &diagram);

  // print one line of segments:
  //
  void print_line(std::vector<int> const &segments_vector);

  std::string diagram_string(time_diagram_t const &diagram);

  // print the arch from a to b with different character.
  //
  void printArch(int a, int b, int k_order, char char1 = '.');
} // namespace inchworm::diagram
