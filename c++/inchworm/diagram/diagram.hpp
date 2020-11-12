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

#include "./../types.hpp"
#include "./../mc/qmc_config_data.hpp"
#include "./utilities.hpp"

namespace inchworm::diagram {

  // Configuration of d(tau) and d_dag(tau') and the split points
  //
  class time_diagram_t {

    // Type containing the operator information and the
    // position in the time ordered list of all d, d_dag
    struct op_t : public fop_t {
      int order_index = 0;
    };

    public:
    std::vector<op_t> op_list;             // list of all operator time ordered
    std::vector<int> split_points;         // position of split points (index of first operator to the right of split time)
    std::vector<fop_t> d_list, d_dag_list; // list of d/d_dag time ordered
    std::vector<int> pos_d;                // position of d in the op_list
    std::vector<int> pos_d_dag;            // idem
    bool is_trivial = true;                // a diagram is considered trivial if no split_times are found between the minimum and maximum tau.

    int perturbation_order() const;
    int size() const;
    double max_tau() const;
    double min_tau() const;

    // Simple function to find the sign of the diagram.
    // Note: this result is the opposite same if we use pos_d_dag
    int sign() const;

    // Constructor
    time_diagram_t(std::vector<fop_t> const &d_list_, std::vector<fop_t> const &d_dag_list_, std::vector<double> const &split_times, int verbose = 0);

    inline time_diagram_t(config_t const &config, std::vector<double> const &split_times, int verbose = 0)
       : time_diagram_t(config.d_list, config.d_dag_list, split_times, verbose) {}
  };
} // namespace inchworm::diagram
