#pragma once

#include "types.hpp"

namespace inchworm {

  // Create an operator containing all possible transitions due to hybridization
  // FIXME: Check actual numerical values of Delta instead??
  many_body_operator create_effective_hyb(gf_struct_t const &gf_struct);

  // Given an atom_diag object and the fundamental operator information, retrieve its block matrix representation
  u_partial_t get_op_block_matrix(atom_diag const &ad, std::string const &bl_name, int idx, bool op_dag);

  // Given the left and right propagator segment insert all operator
  // flavors and take the trace to get the Green function at a given time
  frame_t make_g_frame_from_l_and_r(atom_diag const &ad_imp, gf_struct_t const &gf_struct, u_partial_t const &l, u_partial_t const &r);

  // Initialize bare propagator frame U_0 = exp(-tau H_loc) in the diagonal basis of H_loc
  // CAUTION: Energies are w.r.t. the ground-state energy of the atom_diag object
  frame_t make_bare_u_frame(atom_diag const &ad, double tau);

  // Initialize bare Green function frame
  frame_t make_bare_g_frame(atom_diag const &ad_imp, u_tau_t const &u_tau, gf_struct_t const &gf_struct, double tau_split, double beta);

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

  // Make an exact diagonalization propagator U = Trace_B [exp(-H_bath *(beta-tau)) exp(-H_tot*tau)  ]  /  Trace_B [ exp(-H_bath*beta) ]
  // FIXME Merge with partial_trace_bath
  u_tau_t make_ED_propagator(atom_diag const &ad_tot, atom_diag const &ad_imp, atom_diag const &ad_bath, double beta, int n_tau);

}
