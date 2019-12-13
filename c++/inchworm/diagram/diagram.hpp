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
struct time_and_orbital_t {
  double tau = 0.;
  int orb    = 0;
};

inline bool operator<(time_and_orbital_t const &t1, time_and_orbital_t const &t2) { return (t1.tau < t2.tau); }

// Configuration of c(tau) and cdag(tau') and the split points
//
class time_diagram_t {

  //basic information for an operator
  //
  struct op_t {
    double tau      = 0.;    // time
    bool dag        = false; // true if cdag, false if c
    int orb         = 0;     // orbital index
    int order_index = 0;     // index in the list of all c, cdag
  };

  public:
  std::vector<op_t> op_list;                         // list of all operator time ordered
  std::vector<int> split_points;                     // position of split points
  std::vector<time_and_orbital_t> c_list, cdag_list; // list of c/cdag time ordered
  std::vector<int> pos_c;                            // position of c in the op_list
  std::vector<int> pos_cdag;                         // idem
  bool is_trivial; // a diagram is considered trivial if no split_times are found between the minimum and maximum tau.

  int perturbation_order() const { return c_list.size(); }

  //
  time_diagram_t(std::vector<time_and_orbital_t> const &c, std::vector<time_and_orbital_t> const &cdag, std::vector<double> const &split_times)
     : op_list(2 * c.size()), c_list{c}, cdag_list{cdag} {

    EXPECTS(std::is_sorted(c.begin(), c.end()));
    EXPECTS(std::is_sorted(cdag.begin(), cdag.end()));
    EXPECTS(c.size() == cdag.size());

    split_points.reserve(split_times.size());

    int order = c.size();

    for (int i = 0, j = order; i < order; i++, j++) {
      op_list[i].tau         = c[i].tau;
      op_list[i].orb         = c[i].orb;
      op_list[i].dag         = false;
      op_list[i].order_index = i;

      op_list[j].tau         = cdag[i].tau;
      op_list[j].orb         = cdag[i].orb;
      op_list[j].dag         = true;
      op_list[j].order_index = i;
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
