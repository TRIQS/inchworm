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
//#include <stdio.h>
#include <algorithm>
#include <vector>
#include <numeric>

#include "./utilities.hpp"

// Degrees of freedom of an creation (annihilation) operator:
//
struct time_and_indices_t {
  double tau       = 0.;
  int linear_index = 0; // index as defined in fundamental operators set
};

inline bool operator<(time_and_indices_t const &t1, time_and_indices_t const &t2) { return (t1.tau < t2.tau); }

auto sort_tau = [](auto const &x, auto const &y) { return x.tau < y.tau; };

// Configuration of c(tau) and cdag(tau') and the split points
//
class time_diagram_t {

  //basic information for an operator
  //
  struct op_t {
    double tau       = 0.;    // time
    bool dag         = false; // true if cdag, false if c
    int linear_index = 0;     // index as defined in fundamental operators set
    int order_index  = 0;     // index of the order in time in the list of all c, cdag
  };

  public:
  std::vector<op_t> op_list;                         // list of all operator time ordered
  std::vector<int> split_points;                     // position of split points
  std::vector<time_and_indices_t> c_list, cdag_list; // list of c/cdag time ordered
  std::vector<int> pos_c;                            // position of c in the op_list
  std::vector<int> pos_cdag;                         // idem
  bool is_trivial = true; // a diagram is considered trivial if no split_times are found between the minimum and maximum tau.

  int perturbation_order() const { return c_list.size(); }
  int size() const { return op_list.size(); }
  double max_tau() const { return op_list.back().tau; }
  double min_tau() const { return op_list.front().tau; }

  /*
  bool try_insert_vertices(time_and_indices_t c, time_and_indices_t cdag) {
    // check if different times
    for (int i = 0; i < op_list.size() - 1; i++)
      if ((op_list[i].tau == c.tau) or (op_list[i].tau == cdag.tau)) return false;

    c_list.push_back({c.tau, c.linear_index});
    cdag_list.push_back({cdag.tau, cdag.linear_index});

    std::sort(c_list.begin(), c_list.end(), [](auto const &x, auto const &y) { return x.tau < y.tau; });
    std::sort(cdag_list.begin(), cdag_list.end(), [](auto const &x, auto const &y) { return x.tau < y.tau; });

    reorder();
    return true;
  }

  bool try_insert_split_point() {
     
  }

  time_diagram_t() {}
*/
  //
  time_diagram_t(std::vector<time_and_indices_t> const &c, std::vector<time_and_indices_t> const &cdag, std::vector<double> const &split_times)
     : op_list(2 * c.size()), c_list{c}, cdag_list{cdag} {
//     reorder();
//     }
  
//  bool reorder_op_list() {
    EXPECTS(std::is_sorted(c_list.begin(), c_list.end()));
    EXPECTS(std::is_sorted(cdag_list.begin(), cdag_list.end()));
    EXPECTS(c_list.size() == cdag_list.size());

    split_points.reserve(split_times.size());

    int order = c.size();

    for (int i = 0, j = order; i < order; i++, j++) {
      op_list[i].tau          = c[i].tau;
      op_list[i].linear_index = c[i].linear_index;
      op_list[i].dag          = false;
      op_list[i].order_index  = i;

      op_list[j].tau          = cdag[i].tau;
      op_list[j].linear_index = cdag[i].linear_index;
      op_list[j].dag          = true;
      op_list[j].order_index  = i;
    }

    //if constexpr (verbose) {dd
    //  for (auto op : op_list) { std::printf("%d ", op.order_index); }
    //  std::printf("\n");
    //  for (auto op : op_list) { std::printf("% 4.5f ", op.tau); }
    //  std::printf("\n");
    //}

    std::sort(op_list.begin(), op_list.end(), [](auto const &x, auto const &y) { return x.tau < y.tau; });

    // check that no times are equal (might need to change at some point, rare event, but many Monte Carlo sampling...);
    for (int i = 0; i < op_list.size() - 1; i++) EXPECTS(op_list[i].tau != op_list[i + 1].tau);

    //if constexpr (verbose) {
    //  for (auto op : op_list) { std::printf("%d ", op.order_index); }
    //  std::printf("\n");
    //  for (auto op : op_list) { std::printf("% 4.5f ", op.tau); }
    //  std::printf("\n");
    //}

    if constexpr (verbose) std::printf("split points:\n");

    is_trivial = true; //start by assuming it is trivial and searching for at least one counter example.
    for (auto s_time : split_times) {
      int i = 0;
      for (; i < op_list.size(); i++) {
        EXPECTS(s_time != op_list[i].tau);
        if (s_time < op_list[i].tau) break;
      }
      if (i != 0 and i != op_list.size()) {
        is_trivial = false;
        if constexpr (verbose) std::printf("diagram is not trivial\n");
      }
      split_points.push_back(i);
      if constexpr (verbose) std::printf("%d  % 4.3f\n", i, s_time);
    }

    // posc[i] is the position of the i^th c in op_list (inverse table of order_index)
    for (int i = 0; i < op_list.size(); i++) {
      if (op_list[i].dag)
        pos_cdag.push_back(i);
      else
        pos_c.push_back(i);
    }
  }
};
