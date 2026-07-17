// Copyright (c) 2020--present, The Simons Foundation
// This file is part of inchworm and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./print.hpp"

#include <fmt/core.h>
#include <fmt/color.h>

namespace inchworm::diagram {

  void print_vector(std::vector<int> const &v) {
    for (auto l : v) { std::printf("%d ", l); }
    std::printf("\n");
  }

  static std::vector<fmt::terminal_color> colors{fmt::terminal_color::red,  fmt::terminal_color::yellow, fmt::terminal_color::cyan,
                                                 fmt::terminal_color::blue, fmt::terminal_color::green,  fmt::terminal_color::magenta};

  void print_diag(time_diagram_t const &diagram) {
    for (int j = 0; j < diagram.size(); j++) {
      if (std::any_of(begin(diagram.split_points), end(diagram.split_points), [j](int i) { return i == j + 1; }))
        std::printf(" |");
      else
        std::printf("  ");
    }
    std::printf("\n");
    for (int j = 0; j < diagram.size(); j++) {
      auto &op = diagram.op_list[j];
      if (op.dag)
        fmt::print(fg(colors[op.bl]), "x");
      else
        fmt::print(fg(colors[op.bl]), "o");
      if (j < diagram.size() - 1) std::printf("-");
    }
    std::printf("\n");
  }

  void print_configuration(time_diagram_t const &diagram) {
    print_diag(diagram);
    for (int j = 0; j < diagram.op_list.size(); j++) {
      auto &op = diagram.op_list[j];
      if (op.dag)
        fmt::print("  c†({:.3f})_{}_{}", op.tau, op.bl, op.idx);
      else
        fmt::print("  c({:.3f})_{}_{}", op.tau, op.bl, op.idx);
    }
    std::printf("\n");
  }

  void print_line(std::vector<int> const &segments_vector) {
    int current_segment = 0;
    std::string string1 = "";
    std::string chars   = "  ";
    for (int i = 0; i < segments_vector.size(); i++) {
      if ((current_segment != segments_vector[i]) or (i == 0)) {
        current_segment = segments_vector[i];
        if (current_segment == -1)
          chars = "  ";
        else if (current_segment == 0)
          chars = "--";
        else
          chars = "==";
        if (string1.size() > 0) {
          string1.pop_back();
          string1 += ' ';
        }
      }
      string1 += chars;
    }
    string1.pop_back();
    std::printf("%s", string1.c_str());
  }

  std::string diagram_string(time_diagram_t const &diagram) {
    std::string diagram_order = "";
    for (auto l : diagram.op_list) {
      if (l.dag)
        diagram_order += "x";
      else
        diagram_order += "o";
    }
    return diagram_order;
  }

  void printArch(int a, int b, int k_order, char char1) {
    int ii;
    if (a > b) {
      int tmp = b;
      b       = a;
      a       = tmp;
    } else if (a == b) {
      printf("error a==b");
      exit(1);
    }

    for (ii = 0; ii < 2 * k_order; ii++) {
      if ((ii >= a) and (ii < b))
        printf("%c%c", char1, char1);
      else if (ii == b)
        printf("%c ", char1);
      else
        printf("  ");
    }
    printf("\n");
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
