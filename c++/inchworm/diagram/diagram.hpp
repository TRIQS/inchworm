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

// Function necessary to order different time_and_orbital_t:
//
bool operator<(time_and_orbital_t const &t1, time_and_orbital_t const &t2) { return (t1.tau < t2.tau); }

// Definition of a time_diagram_t:
//
class time_diagram_t {

  struct op_t {
    double tau      = 0.;
    bool dag        = false;
    int orb         = 0;
    int order_index = 0;
  };

  public:
  std::vector<op_t> list;
  std::vector<time_and_orbital_t> c_list, cdag_list;
  //std::vector<double> split_points;
  std::vector<int> split_points;
  std::vector<int> pos_c;
  std::vector<int> pos_cdag;

  int k_order() const { return c_list.size(); }

  time_diagram_t(std::vector<time_and_orbital_t> const &c, std::vector<time_and_orbital_t> const &cdag, std::vector<double> const &split_times)
     : list(2 * c.size()), c_list{c}, cdag_list{cdag} {

    EXPECTS(std::is_sorted(c.begin(), c.end()));
    EXPECTS(std::is_sorted(cdag.begin(), cdag.end()));
    EXPECTS(c.size() == cdag.size());

    split_points.reserve(split_times.size());

    int k_order = c.size();

    for (int i = 0, j = k_order; i < k_order; i++, j++) {
      list[i].tau         = c[i].tau;
      list[i].orb         = c[i].orb;
      list[i].dag         = false;
      list[i].order_index = i;

      list[j].tau         = cdag[i].tau;
      list[j].orb         = cdag[i].orb;
      list[j].dag         = true;
      list[j].order_index = i;
    }

    std::sort(list.begin(), list.end(), [](auto const &x, auto const &y) { return x.tau < y.tau; });

    std::printf("split points:\n");
    for (auto s_time : split_times) {
      int i;
      for (i = 0; i < list.size(); i++) {
        if (s_time < list[i].tau) break;
      }
      split_points.push_back(i);
      std::printf("%d  % 4.3f\n", i, s_time);
    }

    for (int i = 0; i < list.size(); i++) {
      if (list[i].dag)
        pos_cdag.push_back(i);
      else
        pos_c.push_back(i);
    }
  }
};
