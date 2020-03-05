#pragma once
#include "./types.hpp"

#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/utility/macros.hpp>

namespace inchworm {

  using u_frame_t = std::vector<matrix_t>;

  //
  u_frame_t make_zero_propagator_frame(atom_diag const &ad);

  //
  u_frame_t make_bare_propagator_frame(atom_diag const &ad, double tau, bool set_gs_to_0=false);

  // Calculate the Frobenius norm of the u_frame block diagonal matrix:
  double frobenius_norm(u_frame_t const &u_frame);

  // calculate the trace of the u_frame block diagonal matrix:
  double trace(u_frame_t const &u_frame);

  void print(u_frame_t const & u_frame, double factor);
} // namespace inchworm
