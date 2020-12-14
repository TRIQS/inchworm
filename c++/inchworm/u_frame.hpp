#pragma once

#include "types.hpp"

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

  // --------------- General frame / u_partial functionality ---------------

  // Calculate the Frobenius norm of the u_partial block diagonal matrix:
  double frobenius_norm(frame_t const &frame);

  // calculate the trace of the u_frame block diagonal matrix:
  double trace(frame_t const &u_frame);

  // calculate the relative distance of the two frames: ||l - r|| / max(||l||, ||r||)
  double relative_distance(frame_t const &l, frame_t const &r);

  // Create a block_matrix given a gf_struct
  frame_t make_zero_frame(gf_struct_t const &gf_struct);

  // Create a block_matrix given a shape
  // Will be replaced by map(N->matrix{N,N}, shape)
  frame_t make_zero_frame(std::vector<int> const &shape);

  // as_diagonal_block_matrix?
  frame_t make_frame(u_partial_t const &up);

  // as_sparse_block_matrix?
  u_partial_t make_u_partial(frame_t const &u);

  // Multiply two block_matrix
  u_partial_t operator*(u_partial_t const &l, u_partial_t const &r);

  // --------------- Block Gf specific functions -----------------------

  // Get the frame of a Block Green function
  frame_t get_frame(u_tau_t const &u_tau, int idx);

  // Get the frame of a Block Green function using Green function evaluation
  frame_t eval_frame(u_tau_t const &u_tau, double tau);

  // Set a single frame of a Block Green function given a block matrix and the block index
  // TODO Should be as easy as
  //   bgf[bl_][frame_number] << frame[bl_];
  inline void set_frame(frame_t const &frame, u_tau_t &bgf, int frame_number) {
    for (int bl = 0; bl < bgf.size(); bl++) bgf[bl][frame_number] = frame[bl];
  }

  // Calculate the relative distance of two Block Green functions
  // Returns the maximum relative distance of all frames
  double relative_distance(u_tau_t const &l, u_tau_t const &r);

} // namespace inchworm
