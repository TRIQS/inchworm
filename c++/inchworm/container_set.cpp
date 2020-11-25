/*******************************************************************************
 *
 * inchworm: A TRIQS based impurity solver
 *
 * Copyright (c) 2019 The Simons foundation
 *   authors: Nils Wentzell
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
#include "./container_set.hpp"
#include "./u_frame.hpp"
#include "./util.hpp"

namespace inchworm {

  single_step_results_t::single_step_results_t(std::vector<int> const &shape_of_frame) : expansion_order(10, 0), samples_expansion_order(10, 0) {
    frame           = make_zero_frame(shape_of_frame);
    frame_0th_order = frame;
  };

  void single_step_results_t::normalize(double normalization_cte) {
    for (auto &Bl : frame) Bl /= normalization_cte;
    for (auto &Bl : frame_0th_order) Bl /= normalization_cte;
    for (auto &o : expansion_order) o /= normalization_cte;
  };

  void single_step_results_t::print(int verbosity) {

    if (verbosity > 4)
      for (auto Bl : frame) print_matrix(Bl);

    if (verbosity > 3) {
      int max_order = std::min(samples_expansion_order.size(), 15ul);

      std::printf("\n\norder breakdown: \n");
      for (auto k : range(max_order)) std::printf("%10d ", samples_expansion_order[k]);
      std::printf("\n");
      for (auto k : range(max_order)) std::printf("% 10.5f ", expansion_order[k]);
      std::printf("\n");
    }
  }

  void h5_write(h5::group h5group, std::string subgroup_name, container_set const &c) {
    auto grp = h5group.create_group(subgroup_name);
    h5_write(grp, "G_tau", c.G_tau);
  }

  void h5_read(h5::group h5group, std::string subgroup_name, container_set &c) {
    auto grp = h5group.open_group(subgroup_name);
    h5_read(grp, "G_tau", c.G_tau);
  }

} // namespace inchworm
