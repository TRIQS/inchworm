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

namespace inchworm::diagram {
  // Degrees of freedom of an creation (annihilation) operator:
  //
  struct time_and_index_t {
    double tau       = 0.;
    int linear_index = 0; // index as defined in fundamental operators set
  };

  inline bool operator<(time_and_index_t const &t1, time_and_index_t const &t2);

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
    std::vector<op_t> op_list;                       // list of all operator time ordered
    std::vector<int> split_points;                   // position of split points
    std::vector<time_and_index_t> c_list, cdag_list; // list of c/cdag time ordered
    std::vector<int> pos_c;                          // position of c in the op_list
    std::vector<int> pos_cdag;                       // idem
    bool is_trivial = true; // a diagram is considered trivial if no split_times are found between the minimum and maximum tau.

    int perturbation_order() const; 
    int size() const;
    double max_tau() const;
    double min_tau() const;

    // Simple function to find the sign of the diagram. 
    // Note: this result is the opposite same if we use pos_cdag
    int sign() const;

    // Constructor
    time_diagram_t(std::vector<time_and_index_t> const &c, std::vector<time_and_index_t> const &cdag, std::vector<double> const &split_times);
  };
} // namespace inchworm::diagram
