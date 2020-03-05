#pragma once

#include "./u_frame.hpp"
#include "./types.hpp"
#include <triqs/gfs.hpp>
#include <triqs/h5.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/arrays/blas_lapack/dot.hpp>

namespace inchworm {

  /// Return the most significant bit (MSB: the leftmost numbers in a binary representation) via an integer.
  inline uint64_t get_MSB(uint64_t a, int shift) { return (a >> shift); }

  /// Return the least significant bit (LSB: the rigthmost numbers in a binary representation) via an integer.
  inline uint64_t get_LSB(uint64_t a, int shift) { return (a % (1 << shift)); }

  /// Partial sum, tracing over indices above linear_index. Only the linear_index first degrees of freedom will be preserved.
  /**
     * @param ad atom_diag of the system considered here.
     * @param linear_index The linear index (i.e. number) of fundamental operator to be perserved, as defined by the fundamental operator set.
     * @param fct Function to be applied to eigenvalues in atom_diag.
     * @return The partial sum matrix of a function the Hamiltonian.
     */
  triqs::arrays::matrix<double> partial_sum(atom_diag const &ad, int linear_index, std::function<double(double)> fct);
  u_frame_t partial_trace(atom_diag const &ad_full, atom_diag const &ad_target, std::function<double(double)> fct);
  scalar_t trace(atom_diag const &ad_full, std::function<double(double)> fct);

  void print_energies(std::vector<std::vector<double>> const &E);
  void print_eigensystems(atom_diag const &ad); 
  void print_matrix(triqs::arrays::matrix<double> m);

} // namespace inchworm
