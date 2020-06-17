#pragma once

#include <numeric>
#include <bitset>
#include <iomanip>

#include "./u_frame.hpp"
#include "./types.hpp"
#include <triqs/gfs.hpp>
#include <triqs/h5.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/arrays/blas_lapack/dot.hpp>

namespace inchworm {

  void print_binary(unsigned int n, int total_bits);

  /// Return the most significant bit (MSB: the leftmost numbers in a binary representation) via an integer.
  uint64_t get_MSB(uint64_t a, int shift);

  /// Return the least significant bit (LSB: the rigthmost numbers in a binary representation) via an integer.
  uint64_t get_LSB(uint64_t a, int shift);
 
  // function that calls partial_trace_bath. FIXME Could be merge at some point, probably no need to have two seperates functions.
  //u_tau_t make_ED_propagator(atom_diag const &ad_tot, atom_diag const &ad_atom, atom_diag const &ad_bath, double beta, int n_tau);

  /// Partial trace, tracing over indices above linear_index, i.e. the bath degree of freedom. Only the linear_index first degrees of freedom will be preserved.
  /**
     * @param ad_tot atom_diag of the full Hamiltonian H = H_loc + H_bath + H_hyb.
     * @param ad_imp  atom_diag of the local Hamiltonian H_loc.
     * @param ad_bath atom_diag of the bath Hamiltonian H_bath (in the full basis).
     * @param beta inverse temperature.
     * @param tau 0 < tau < beta.
     * @return The partial sum matrix of a function H. The result is a block diagonal matrix, with blocks and indices in the same order as in ad_imp.
     */
  u_frame_t partial_trace_bath(atom_diag const &ad_tot, atom_diag const &ad_target, atom_diag const &ad_bath, double beta, double dtau);

  // Trace over all degrees of freedom of atom_diag.
  scalar_t trace(atom_diag const &ad_tot, std::function<double(double)> fct);

  // printing:
  void print_energies(std::vector<std::vector<double>> const &E);
  void print_eigensystems(atom_diag const &ad);
  void print_atom_diag(atom_diag const &ad);
  void print_matrix(triqs::arrays::matrix<double> m, scalar_t factor = 1.0);
  void fprint(u_tau_t u_tau, int N_tau);
  void print(u_tau_t u_tau, int frame_number);
  void print(u_tau_t u_tau, double tau);
  void print(u_frame_t u_frame);
  void print(u_partial_t u_partial);

  // necessary to fill the propagator at each step of the inchworm:
  template<typename T>
  void assign_u_frame_to_propagator(block_gf<imtime, T> &u_tau, u_frame_t const &u_frame, int frame_number, scalar_t factor = 1.0) {
    for (int bl = 0; bl < u_tau.size(); bl++) u_tau[bl][frame_number] = factor * u_frame[bl];
    return;
  }


} // namespace inchworm
