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
#include "./diagram.hpp"

namespace inchworm::diagram {

  int time_diagram_t::perturbation_order() const { return d_list.size(); }

  int time_diagram_t::size() const { return op_list.size(); }

  double time_diagram_t::min_tau() const { return op_list.front().tau; }

  double time_diagram_t::max_tau() const { return op_list.back().tau; }

  int time_diagram_t::sign() const {
    // Note: use pos_d_dag, not pos_d as we consider <d_k ddag_k .. d_0 ddag_0> from the right
    return (std::accumulate(pos_d_dag.begin(), pos_d_dag.end(), 0) % 2 == 0 ? 1 : -1);
  }

  time_diagram_t::time_diagram_t(std::vector<fop_t> const &d_list_, std::vector<fop_t> const &d_dag_list_, std::vector<double> const &split_times,
                                 int verbose)
     : op_list(2 * d_list_.size()), d_list{d_list_}, d_dag_list{d_dag_list_} {
    EXPECTS(d_list.size() == d_dag_list.size());

    int order = d_list.size();
    if (order == 0) return;

    // Initialize the time-ordered list of all operators
    for (int i = 0, j = order; i < order; i++, j++) {
      op_list[i].tau          = d_list[i].tau;
      op_list[i].dag          = false;
      op_list[i].linear_index = d_list[i].linear_index;
      op_list[i].order_index  = i;

      op_list[j].tau          = d_dag_list[i].tau;
      op_list[j].dag          = true;
      op_list[j].linear_index = d_dag_list[i].linear_index;
      op_list[j].order_index  = i;
    }
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
