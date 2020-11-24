#pragma once

#include <numeric>
#include <bitset>
#include <iomanip>

#include "./u_frame.hpp"
#include "./types.hpp"
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <h5/h5.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/atom_diag/atom_diag.hpp>
#include <nda/blas/dot.hpp>

namespace inchworm {

  // printing:
  void print_energies(std::vector<std::vector<double>> const &E);
  void print_eigensystems(atom_diag const &ad);
  void print_atom_diag(atom_diag const &ad);
  void print_matrix(triqs::arrays::matrix<double> m, scalar_t factor = 1.0);
  void fprint(u_tau_t u_tau, int N_tau);
  void print(u_tau_t u_tau, int frame_number);
  void print(u_tau_t u_tau, double tau);
  void print(frame_t u_frame);
  void print(u_partial_t u_partial);

  //------------------------------

  // Create an operator containing all possible transitions due to hybridization
  // FIXME: Check actual numerical values of Delta instead??
  many_body_operator create_effective_hyb(gf_struct_t const &gf_struct);

  /// =============== TODO: Move into testing functionality =================

  /**
   * Calculate U(tau) = Tr_bath [ exp[-(beta - tau) * H_bath] * exp(-tau * H_tot) ] / Z_bath
   *
   * Within atom_diag, perform partial trace over indices above nfops_imp, i.e. the bath degrees of freedom.
   * Only the nfops_imp first degrees of freedom will be preserved.
   *
   * @param ad_tot atom_diag of the full Hamiltonian H = H_loc + H_bath + H_hyb.
   * @param ad_imp  atom_diag of the local Hamiltonian H_loc.
   * @param ad_bath atom_diag of the bath Hamiltonian H_bath (in the full basis).
   * @param beta inverse temperature.
   * @param tau 0 < tau < beta.
   * @return The partial sum matrix of a function H. The result is a block diagonal matrix, with blocks and indices in the same order as in ad_imp.
   */
  frame_t partial_trace_bath(atom_diag const &ad_tot, atom_diag const &ad_target, atom_diag const &ad_bath, double beta, double tau);

  /// =============== End: Move into testing functionality =================

} // namespace inchworm
