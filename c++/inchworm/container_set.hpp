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
#pragma once

#include "types.hpp"

namespace inchworm {

  // structure to gather result of one Monte Carlo run:
  struct qmc_results_t {

    frame_t frame;
    frame_t frame_0th_order;

    std::vector<frame_t> frame_by_order;
    std::vector<double> order_histogram;

    double average_order = 0.0;

    int status = 0;

    double auto_corr_time;

    qmc_results_t(std::vector<int> const &shape_of_frame);

    void normalize(double normalization_cte);
    void print(int verbosity = 4);
  };

  /// The collection of all output containers in solver_core
  struct container_set {

    /// Greens function in imaginary time
    g_tau_t G_tau;

    /// Order-resolved propagator
    std::vector<g_tau_t> G_tau_by_order;

    /// Order-resolved propagator
    std::vector<u_tau_t> u_tau_by_order;

    /// Order histograms
    std::vector<std::vector<double>> order_histograms;

    /// Function that writes all containers to hdf5 file
    friend void h5_write(h5::group h5group, std::string subgroup_name, container_set const &c);

    /// Function that reads all containers from hdf5 file
    friend void h5_read(h5::group h5group, std::string subgroup_name, container_set &c);
  };

} // namespace inchworm
