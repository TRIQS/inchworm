// Copyright (c) 2019--present, The Simons Foundation
// This file is part of inchworm and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once

#include "types.hpp"

namespace inchworm {

  // structure to gather result of one Monte Carlo run:
  struct qmc_results_t {

    frame_t frame;
    scalar_t weight_zeroth_order = 0.0;
    std::vector<scalar_t> errs_frame;

    std::vector<frame_t> frame_by_order;
    std::vector<scalar_t> err_frame_by_order;

    std::vector<double> order_histogram;

    scalar_t average_sign;
    double average_order = 0.0;

    int status = 0;

    double auto_corr_time;

    qmc_results_t(std::vector<int> const &shape_of_frame);

    void normalize(scalar_t normalization_cte);

    void print(int verbosity = 4);
  };

  /// The collection of all output containers in solver_core
  struct container_set {

    /// Greens function in imaginary time
    g_tau_t G_tau;

    /// Last obtained frame error block-resolved
    std::vector<scalar_t> errs_frame;

    /// Order-resolved propagator
    std::vector<g_tau_t> G_tau_by_order;

    /// Order-resolved propagator
    std::vector<u_tau_t> u_tau_by_order;

    /// Last obtained frame error order-resolved
    std::vector<scalar_t> err_frame_by_order;

    /// Order histograms
    std::vector<std::vector<double>> order_histograms;

    /// Function that writes all containers to hdf5 file
    friend void h5_write(h5::group h5group, std::string subgroup_name, container_set const &c);

    /// Function that reads all containers from hdf5 file
    friend void h5_read(h5::group h5group, std::string subgroup_name, container_set &c);
  };

} // namespace inchworm
