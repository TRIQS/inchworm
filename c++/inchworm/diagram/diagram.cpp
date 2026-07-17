// Copyright (c) 2020--present, The Simons Foundation
// This file is part of inchworm and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./diagram.hpp"
#include "./proper_enum.hpp"

namespace inchworm::diagram {

  int time_diagram_t::perturbation_order() const { return d_list.size(); }

  int time_diagram_t::size() const { return op_list.size(); }

  double time_diagram_t::min_tau() const { return op_list.front().tau; }

  double time_diagram_t::max_tau() const { return op_list.back().tau; }

  int time_diagram_t::sign() const { return _sign; }

  time_diagram_t::time_diagram_t(std::vector<fop_t> const &d_list_, std::vector<fop_t> const &d_dag_list_, std::vector<double> const &split_times,
                                 int verbose)
     : op_list(2 * d_list_.size()), d_list{d_list_}, d_dag_list{d_dag_list_} {
    EXPECTS(d_list.size() == d_dag_list.size());

    int order = d_list.size();
    if (order == 0) return;

    // Initialize the time-ordered list of all operators
    for (int i = 0; i < order; i++) {
      op_list[2 * i].tau          = d_dag_list[i].tau;
      op_list[2 * i].dag          = true;
      op_list[2 * i].linear_index = d_dag_list[i].linear_index;
      op_list[2 * i].bl           = d_dag_list[i].bl;
      op_list[2 * i].idx          = d_dag_list[i].idx;
      op_list[2 * i].order_index  = i;

      op_list[2 * i + 1].tau          = d_list[i].tau;
      op_list[2 * i + 1].dag          = false;
      op_list[2 * i + 1].linear_index = d_list[i].linear_index;
      op_list[2 * i + 1].bl           = d_list[i].bl;
      op_list[2 * i + 1].idx          = d_list[i].idx;
      op_list[2 * i + 1].order_index  = i;
    }

    std::vector<int> ivec(op_list.size());
    std::iota(begin(ivec), end(ivec), 0);
    std::sort(begin(ivec), end(ivec), [this](int i, int j) { return op_list[i] < op_list[j]; });
    _sign = find_parity(ivec);

    std::sort(op_list.begin(), op_list.end(), std::less<>{});
    for (int i = 0; i < op_list.size() - 1; i++) EXPECTS(op_list[i].tau != op_list[i + 1].tau);

    // For each split_time determine the op_list index of the operator to the right (i.e. the split point)
    if (verbose > 3) std::printf("split points:\n");
    split_points.reserve(split_times.size());
    for (auto s_time : split_times) {
      int i = 0;
      for (; i < op_list.size(); i++) {
        EXPECTS(s_time != op_list[i].tau);
        if (s_time < op_list[i].tau) break;
      }
      split_points.push_back(i);
      if (verbose > 3) std::printf("%d  % 4.3f\n", i, s_time);
    }

    // Diagram is trivial if no split_times are found between the smallest and largest operator time
    is_trivial = std::all_of(cbegin(split_points), cend(split_points), [&](int i) { return i == 0 or i == op_list.size(); });
    if (is_trivial and verbose > 3) std::printf("diagram is trivial\n");

    // For each d and d_dag, find its position in op_list (inverse table of order_index)
    for (int i = 0; i < op_list.size(); i++) {
      if (op_list[i].dag)
        pos_d_dag.push_back(i);
      else
        pos_d.push_back(i);
    }
  }
} // namespace inchworm::diagram
