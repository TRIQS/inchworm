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

  inline bool operator<(time_and_index_t const &t1, time_and_index_t const &t2) { return (t1.tau < t2.tau); }

  //auto sort_tau = [](auto const &x, auto const &y) { return x.tau < y.tau; };

  int time_diagram_t::perturbation_order() const { return d_list.size(); }
  int time_diagram_t::size() const { return op_list.size(); }
  double time_diagram_t::max_tau() const { return op_list.back().tau; }
  double time_diagram_t::min_tau() const { return op_list.front().tau; }

  // Simple function to find the sign of the diagram.
  // Note: important to use pos_d_dag, not pos_d
  int time_diagram_t::sign() const { return (std::accumulate(pos_d_dag.begin(), pos_d_dag.end(), 0) % 2 == 0 ? 1 : -1); }

  //
  time_diagram_t::time_diagram_t(std::vector<time_and_index_t> const &d, std::vector<time_and_index_t> const &d_dag,
                                 std::vector<double> const &split_times, int verbose)
     : op_list(2 * d.size()), d_list{d}, d_dag_list{d_dag} {

    std::sort(d_dag_list.begin(), d_dag_list.end(), [](auto const &x, auto const &y) { return x.tau < y.tau; });
    std::sort(d_list.begin(), d_list.end(), [](auto const &x, auto const &y) { return x.tau < y.tau; });

    EXPECTS(d_list.size() == d_dag_list.size());
    split_points.reserve(split_times.size());

    int order = d_list.size();

    if (d.size() == 0) {
      //std::printf("warning: zero lenght! \n");
      //fflush(stdout);
      return;
    }
    for (int i = 0, j = order; i < order; i++, j++) {
      op_list[i].tau          = d_list[i].tau;
      op_list[i].linear_index = d_list[i].linear_index;
      op_list[i].dag          = false;
      op_list[i].order_index  = i;

      op_list[j].tau          = d_dag_list[i].tau;
      op_list[j].linear_index = d_dag_list[i].linear_index;
      op_list[j].dag          = true;
      op_list[j].order_index  = i;
    }

    std::sort(op_list.begin(), op_list.end(), [](auto const &x, auto const &y) { return x.tau < y.tau; });

    // check that no times are equal (might need to change at some point, rare event, but many Monte Carlo sampling...);
    for (int i = 0; i < op_list.size() - 1; i++) EXPECTS(op_list[i].tau != op_list[i + 1].tau);

    if(verbose > 3) std::printf("split points:\n");

    is_trivial = true; //start by assuming it is trivial and searching for at least one counter example.
    for (auto s_time : split_times) {
      int i = 0;
      for (; i < op_list.size(); i++) {
        EXPECTS(s_time != op_list[i].tau);
        if (s_time < op_list[i].tau) break;
      }
      if (i != 0 and i != op_list.size()) {
        is_trivial = false;
        if (verbose > 3) std::printf("diagram is not trivial\n");
      }
      split_points.push_back(i);
      if (verbose > 3) std::printf("%d  % 4.3f\n", i, s_time);
    }
    // posc[i] is the position of the i^th c in op_list (inverse table of order_index)
    for (int i = 0; i < op_list.size(); i++) {
      if (op_list[i].dag)
        pos_d_dag.push_back(i);
      else
        pos_d.push_back(i);
    }
  }
} // namespace inchworm::diagram
