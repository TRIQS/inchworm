// Copyright (c) 2019--present, The Simons Foundation
// This file is part of inchworm and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

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
