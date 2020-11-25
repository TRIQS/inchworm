#pragma once

#include "types.hpp"

namespace inchworm {

  void print_energies(std::vector<std::vector<double>> const &E);
  void print_eigensystems(atom_diag const &ad);
  void print_atom_diag(atom_diag const &ad);
  void print_matrix(triqs::arrays::matrix<double> m, scalar_t factor = 1.0);
  void fprint(u_tau_t u_tau, int N_tau);
  void print(u_tau_t u_tau, int frame_number);
  void print(u_tau_t u_tau, double tau);
  void print(frame_t u_frame);
  void print(u_partial_t u_partial);

} // namespace inchworm
