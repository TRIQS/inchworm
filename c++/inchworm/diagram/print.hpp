// Copyright (c) 2019--present, The Simons Foundation
// This file is part of inchworm and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once

#include "./diagram.hpp"
#include "./segment.hpp"

namespace inchworm::diagram {

  void print_vector(std::vector<int> const &v);

  // Print diagram and its split point above.
  void print_diag(time_diagram_t const &diagram);

  void print_configuration(time_diagram_t const &diagram);

  // print one line of segments:
  void print_line(std::vector<int> const &segments_vector);

  std::string diagram_string(time_diagram_t const &diagram);

  // print the arch from a to b with different character.
  void printArch(int a, int b, int k_order, char char1 = '.');

  // Print a single segment
  void print_segment(segment_t const &segment, time_diagram_t const &diagram);

  // Print a set of segments
  void print_set(std::vector<segment_t> const &segment_list, set_of_segments_t const &set, time_diagram_t const &diagram);

} // namespace inchworm::diagram
