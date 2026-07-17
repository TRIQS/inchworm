// Copyright (c) 2019--present, The Simons Foundation
// This file is part of inchworm and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once

#include "../types.hpp"
#include "../config.hpp"

#include <vector>

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
    std::vector<op_t> op_list;                    // Time ordered list of all operators
    std::vector<int> split_points;                // Position of split points (index of first operator to the right of split time)
    std::vector<fop_t> const &d_list, d_dag_list; // Time ordered list of d/d_dag
    bool is_trivial = true;                       // Diagram is trivial if no split_times are found between the smallest and largest operator time
    std::vector<int> pos_d;                       // Position of d in the op_list
    std::vector<int> pos_d_dag;                   // Position of d_dag in the op_list

    /// The perturbation order
    int perturbation_order() const;

    /// The total number of operators in the diagram
    int size() const;

    /// The smallest operator time
    double min_tau() const;

    /// The largest operator time
    double max_tau() const;

    // The sign of the diagram arising from the order of operators
    int sign() const;

    // Constructor
    time_diagram_t(std::vector<fop_t> const &d_list_, std::vector<fop_t> const &d_dag_list_, std::vector<double> const &split_times, int verbose = 0);

    time_diagram_t(config_t const &config, std::vector<double> const &split_times, int verbose = 0)
       : time_diagram_t(config.d_list, config.d_dag_list, split_times, verbose) {}

    private:
    int _sign = 1; // Order zero default
  };

} // namespace inchworm::diagram
