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

  qmc_results_t::qmc_results_t(std::vector<int> const &shape_of_frame)
     : frame{make_zero_frame(shape_of_frame)}, errs_frame(shape_of_frame.size(), 0.0) {};

  void qmc_results_t::normalize(scalar_t normalization_cte) {
    if (normalization_cte == scalar_t{0.}) TRIQS_RUNTIME_ERROR << "Error in result normalization: Division by zero.";
    frame /= normalization_cte;
    weight_zeroth_order /= normalization_cte;
    for (auto &err_frame_bl : errs_frame) err_frame_bl /= normalization_cte;
    for (auto &frame_k : frame_by_order) frame_k /= normalization_cte;
  };

  void qmc_results_t::print(int verbosity) {

    if (verbosity > 0) {
      std::printf("     max(abs(frame)): %.4e\n", max_element(nda::map([](matrix_t const &m) { return max_element(abs(m)); })(frame)));
      std::printf("     min(abs(frame)): %.4e\n", min_element(nda::map([](matrix_t const &m) { return min_element(abs(m)); })(frame)));
      std::printf("     average_sign: %4f\n", average_sign);
      std::printf("     average_order: %4f\n", average_order);
      if (order_histogram.size() > 0) {
        std::printf("     order_histogram: [");
        for (auto v : order_histogram) std::printf(" %.3f, ", v);
        std::printf("]\n");
      }
      if (frame_by_order.size() > 0) {
        std::printf("     order contribution: [");
        for (auto &frame_k : frame_by_order) std::printf(" %.3f, ", norm(frame_k) / norm(frame));
        std::printf("]\n");
      }
    }

    if (verbosity > 4)
      for (auto Bl : frame) print_matrix(Bl);
  }

  void h5_write(h5::group h5group, std::string subgroup_name, container_set const &c) {
    auto grp = h5group.create_group(subgroup_name);
    h5_write(grp, "G_tau", c.G_tau);
    h5_write(grp, "errs_frame", c.errs_frame);
    h5_write(grp, "G_tau_by_order", c.G_tau_by_order);
    h5_write(grp, "u_tau_by_order", c.u_tau_by_order);
    h5_write(grp, "err_frame_by_order", c.err_frame_by_order);
    h5_write(grp, "order_histograms", c.order_histograms);
  }

  void h5_read(h5::group h5group, std::string subgroup_name, container_set &c) {
    auto grp = h5group.open_group(subgroup_name);
    h5_read(grp, "G_tau", c.G_tau);
    h5::try_read(grp, "errs_frame", c.errs_frame);
    h5::try_read(grp, "G_tau_by_order", c.G_tau_by_order);
    h5::try_read(grp, "u_tau_by_order", c.u_tau_by_order);
    h5::try_read(grp, "err_frame_by_order", c.err_frame_by_order);
    h5_read(grp, "order_histograms", c.order_histograms);
  }

} // namespace inchworm
