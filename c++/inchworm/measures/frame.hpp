// Copyright (c) 2019--present, The Simons Foundation
// This file is part of inchworm and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once

#include "../params.hpp"
#include "../config.hpp"
#include "../container_set.hpp" // qmc_results_t

namespace inchworm::measures {

  /// Measurement of a single frame of a propagator or Green function
  struct frame {

    frame(params_t const &, config_t const &config, frame_t const &frame, qmc_results_t &results);

    /// Invoke a single measurement
    void accumulate(scalar_t);

    /// Finalize
    void collect_results(mpi::communicator const &comm);

    private:
    // Print verbosity
    int verbosity;

    // The Monte-Carlo configuration
    config_t const &config;

    // The current frame
    frame_t const &curr_frame;

    // References to the accumulation frames
    frame_t &acc_frame;
    scalar_t &weight_zeroth_order;

    // Errors of the (0,0) component for each frame[bl]
    std::vector<scalar_t> &errs_frame;

    // The scalar accumulator for the auto-correlation analysis
    log_binning<scalar_t> log_acc = {0.0, -1};

    // The scalar accumulator for the error analysis
    std::vector<lin_binning<scalar_t>> lin_acc;

    // The number of samples
    long long N_samples = 0;
  };

} // namespace inchworm::measures
