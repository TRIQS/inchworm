// Copyright (c) 2020--present, The Simons Foundation
// This file is part of inchworm and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./average_order.hpp"

namespace inchworm::measures {

  average_order::average_order(params_t const &, config_t const &config, qmc_results_t &results)
     : config(config), average_order_ref(results.average_order) {
    average_order_ref = 0.0;
  }

  void average_order::accumulate(scalar_t) {
    average_order_ref += config.size();
    ++N;
  }

  void average_order::collect_results(mpi::communicator const &comm) {
    N = mpi::all_reduce(N, comm);

    // Reduce and normalize
    average_order_ref = mpi::all_reduce(average_order_ref, comm);
    average_order_ref = average_order_ref / N;
  }

} // namespace inchworm::measures
