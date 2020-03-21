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
#include "./print.hpp"

///////////////////////////////////
////////// PRINT //////////////////
///////////////////////////////////

namespace inchworm::diagram {
  void print_vector(std::vector<int> const &v) {
    for (auto l : v) { std::printf("%d ", l); }
    std::printf("\n");
  }

  // print diagram and its split point above.
  //
  void print_diag(time_diagram_t const &diagram) {
    for (int j = 0; j < diagram.op_list.size(); j++) {
      if (std::any_of(begin(diagram.split_points), end(diagram.split_points), [j](int i) { return i == j + 1; }))
        std::printf(" |");
      else
        std::printf("  ");
    }
    std::printf("\n");
    for (int j = 0; j < diagram.op_list.size(); j++) {
      if (diagram.op_list[j].dag)
        std::printf("x");
      else
        std::printf("o");
      if (j < diagram.op_list.size() - 1) std::printf("-");
    }
    std::printf("\n");
  }

  void print_configuration(time_diagram_t const &diagram) {
    print_diag(diagram);
    for (int j = 0; j < diagram.op_list.size(); j++) {
      if (diagram.op_list[j].dag)
        std::printf(" c(%2.5f)_%d", diagram.op_list[j].tau, diagram.op_list[j].linear_index);
      else
        std::printf(" c^+(%2.5f)_%d", diagram.op_list[j].tau, diagram.op_list[j].linear_index);
    }
    std::printf("\n");
  }

  // print one line of segments:
  //
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

  // print the arch from a to b with different character.
  //
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
      else if ((ii == b))
        printf("%c ", char1);
      else
        printf("  ");
    }
    printf("\n");
  }
} // namespace inchworm::diagram
