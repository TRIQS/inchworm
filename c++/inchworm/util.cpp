// Copyright (c) 2020--present, The Simons Foundation
// This file is part of inchworm and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./util.hpp"
#include <iomanip>

namespace inchworm {

  void print_energies(std::vector<std::vector<double>> const &E) {
    for (auto sp : E) {
      for (auto l : sp) { std::printf("% 2.3f ", l); }
      std::printf("\n");
    }
    std::printf("\n");
  }

  void print_eigensystems(atom_diag const &ad) {
    for (auto sp : ad.get_eigensystems()) {
      for (auto l : sp.eigenvalues) { std::printf("% 2.3f ", l); }
      std::printf("\n\n");

      for (int i = 0; i < sp.eigenvalues.size(); i++) {
        for (int j = 0; j < sp.eigenvalues.size(); j++) { std::printf("% 2.3f ", sp.unitary_matrix(i, j)); }
        std::printf("\n");
      }
      //for (auto u : sp.unitary_matrix) { TRIQS_PRINT(u); }
      std::printf("\n\n");
    }
    std::printf("\n");
  }

  void print_binary(unsigned int n, int total_bits) {

    for (int i = 0; i < total_bits; i++) {
      int bit = (1 << (total_bits - i - 1));
      if ((n & bit) != 0)
        std::printf("1");
      else
        std::printf("0");
    }
    std::printf(" ");
  }

  void print_atom_diag(atom_diag const &ad) {

    int N = 0;
    for (int s = 0; s < ad.n_subspaces(); s++) {
      auto sp     = ad.get_eigensystems()[s];
      auto E_Udag = matrix_t{dagger(sp.unitary_matrix)};
      for (int i = 0; i < sp.eigenvalues.size(); i++) {
        for (int j = 0; j < sp.eigenvalues.size(); j++) E_Udag(i, j) *= sp.eigenvalues[i];
      }
      auto H = sp.unitary_matrix * E_Udag;

      std::printf("\nblock: %d \n   ", N++);
      auto fss = ad.get_fock_states()[s];
      for (int i = 0; i < ad.get_subspace_dim(s); i++) {
        printf("  ");
        print_binary(fss[i], ad.get_fops().data().size());
      }
      std::printf("\n");
      print_matrix(H);
    }
    std::printf("\n");
  }

  void print_matrix(matrix_t const &m, double factor) {

    for (int i = 0; i < m.shape()[0]; i++) {
      std::printf("\n [");
      for (int j = 0; j < m.shape()[1]; j++) {
        if (std::abs(m(i, j)) == 0.0)
          std::printf("  .        ");
        else
          std::printf("% 10.3e ", m(i, j) * factor);
      }
      std::printf("]");
    }
    std::printf("\n");
  }

  void fprint(u_tau_t const &u_tau, int N_tau) {
    FILE *f = fopen("u_tau.dat", "w");
    for (int i_tau = 0; i_tau < N_tau; i_tau++) {
      fprintf(f, "%d  ", i_tau);
      for (int bl = 0; bl < u_tau.size(); bl++) print_matrix((matrix_t)u_tau[bl][i_tau]);
      fprintf(f, "\n");
    }
    fclose(f);
    return;
  }

  void print(u_partial_t const &u_partial) {
    for (auto &[bl, mat] : u_partial) {
      std::printf("\n\nblock: %d\n", bl);
      print_matrix(mat);
    }
    std::cout << "\n";
    return;
  }

  void print(frame_t const &u_frame) {
    for (int bl = 0; bl < u_frame.size(); bl++) {
      std::printf("\n\nblock: %d\n", bl);
      print_matrix(u_frame[bl]);
    }
    std::cout << "\n";
    return;
  }

  void print(u_tau_t const &u_tau, double tau) {
    for (int bl = 0; bl < u_tau.size(); bl++) { print_matrix((matrix_t)u_tau[bl](tau)); }
    std::cout << "\n";
    return;
  }

  void print(u_tau_t const &u_tau, int frame_number) {
    for (int bl = 0; bl < u_tau.size(); bl++) { std::cout << std::setprecision(10) << u_tau[bl][frame_number]; }
    std::cout << "\n";
    return;
  }

} // namespace inchworm
