// Copyright (c) 2020--present, The Simons Foundation
// This file is part of inchworm and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once

#include "types.hpp"

namespace inchworm {

  void print_energies(std::vector<std::vector<double>> const &E);
  void print_eigensystems(atom_diag const &ad);
  void print_atom_diag(atom_diag const &ad);
  void print_matrix(matrix_t const &m, double factor = 1.0);
  void fprint(u_tau_t const &u_tau, int N_tau);
  void print(u_tau_t const &u_tau, int frame_number);
  void print(u_tau_t const &u_tau, double tau);
  void print(frame_t const &u_frame);
  void print(u_partial_t const &u_partial);

} // namespace inchworm
