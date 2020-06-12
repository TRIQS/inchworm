#pragma once
#include "./types.hpp"

#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/utility/macros.hpp>
#include <itertools/itertools.hpp>

#include <algorithm>

namespace inchworm {

  using u_frame_t = std::vector<matrix_t>;

  using u_partial_t = std::vector<std::pair<int, matrix_t>>;

  //
  u_frame_t make_zero_propagator_frame(atom_diag const &ad);

  //
  u_frame_t make_bare_propagator_frame(atom_diag const &ad, double tau, bool set_gs_to_0 = false);

  u_frame_t make_u_frame(u_partial_t const &up);

  u_partial_t make_u_partial(u_frame_t const &u);

  void multiply_assign(u_frame_t &u_frame, u_partial_t const &l, u_partial_t const &r);
  //u_frame_t operator*(u_partial_t const &l, u_partial_t const &r);

  // Calculate the Frobenius norm of the u_frame block diagonal matrix:
  double frobenius_norm(u_frame_t const &u_frame);

  // calculate the trace of the u_frame block diagonal matrix:
  double trace(u_frame_t const &u_frame);

  void print(u_frame_t const &u_frame, double factor);
} // namespace inchworm
