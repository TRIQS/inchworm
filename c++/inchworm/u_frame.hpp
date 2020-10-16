#pragma once
#include "./types.hpp"

#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/utility/macros.hpp>
#include <itertools/itertools.hpp>

#include <algorithm>

namespace inchworm {

  // TODO Change type of u_partial_t to
  // vector<pair<int, matrix_t>>;
  // -> array<pair<int, matrix>> ???
  // -> struct{ array<int>, array<matrix> }; ???
  //
  // Rename to:
  // - sparse_block_matrix_t ??
  // 
  // Move it into atom diag??
  // Atom diag should have an internal operator repr.
  // which allows for op->sparse_block_matrix
  // Should this apply for just c / c_dag??
  // Or also also generic operators ??  
  // (cdag_i c_j + cdag_k .. )-> full_block_matrix ???

  // TODO Change type of frame_t to
  // -> array<matrix_t>;
  // 
  // Rename frame_t to:
  // - diag_block_matrix_t ??

  // as_diagonal_block_matrix
  frame_t make_u_frame(u_partial_t const &up);

  // as_sparse_block_matrix
  u_partial_t make_u_partial(frame_t const &u);

  // Create an empty frame (block diagonal matrix: vector of matrix_t)
  frame_t make_zero_propagator_frame(atom_diag const &ad);

  // initialize bare propagator frame U_0 = exp(-tau H_loc) in the diagonal basis of H_loc
  frame_t make_bare_propagator_frame(atom_diag const &ad, double tau, bool set_gs_to_0 = false);

  frame_t make_bare_g_frame(atom_diag const &ad_imp, u_tau_t const &u_tau, std::map<int, std::pair<int, int>> const &map_lin_idx_to_block_inner,
                            gf_struct_t const &gf_struct, double tau_split, double beta);

  // Create a diag_block_matrix given a gf_struct
  frame_t make_frame(gf_struct_t const & gf_struct);

  // Create a diag_block_matrix given a shape
  // Will be replaced by map(N->matrix{N,N}, shape)
  frame_t make_frame(std::vector<long> const & shape);

  // TODO Keep it with the u_partial_t implementation
  u_partial_t operator*(u_partial_t const &l, u_partial_t const &r);

  // Apply, to a non-diagonal propagator, a c or cdag operator from the right
  u_partial_t apply_op_from_right(u_partial_t const &l, int lin_index, bool op_dag, atom_diag const & ad);

  // Calculate the Frobenius norm of the u_frame block diagonal matrix:
  double frobenius_norm(u_partial_t const &u_partial);

  // Calculate the Frobenius norm of the u_partial block diagonal matrix:
  double frobenius_norm(frame_t const &frame);

  // calculate the trace of the u_frame block diagonal matrix:
  double trace(frame_t const &u_frame);

  // calculate the trace of the u_partial block non-diagonal matrix:
  inline double trace(u_partial_t const &u_partial) { return trace(make_u_frame(u_partial)); }

  // Get the frame of a Green function
  // Only used in one place -> remove
  frame_t get_frame(u_tau_t const & u_tau, int idx);

  // calculate the relative distance of the two frames: ||l - r|| / max(||l||, ||r||)
  double relative_distance(frame_t const &l, frame_t const &r);

  // calculate the relative distance of two propagators
  // maximum relative distance of all frames
  double relative_distance(u_tau_t const &l, u_tau_t const &r);

  void print(frame_t const &u_frame, double factor);

  // Given the left and right propagator segment insert all operator 
  // flavors and take the trace to get the Green function at a given 
  // time bl == -1 is translated into zero-initialized matrices
  frame_t make_g_frame_from_l_and_r(atom_diag const &ad_imp, std::map<int, std::pair<int, int>> const &map_lin_idx_to_block_inner,
                                    gf_struct_t const &gf_struct, u_partial_t const &l, u_partial_t const &r);

} // namespace inchworm
