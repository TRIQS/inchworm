#pragma once
#include "./types.hpp"

#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/utility/macros.hpp>
#include <itertools/itertools.hpp>

#include <algorithm>

namespace inchworm {

  using frame_t = std::vector<matrix_t>;

  using u_partial_t = std::vector<std::pair<int, matrix_t>>;

  using g_frame_t = std::vector<matrix_t>;

  //
  frame_t make_zero_propagator_frame(atom_diag const &ad);

  //
  frame_t make_bare_propagator_frame(atom_diag const &ad, double tau, bool set_gs_to_0 = false);

  frame_t make_bare_g_frame(atom_diag const &ad_imp, u_tau_t const &u_tau, std::map<int, std::pair<int, int>> const &map_lin_idx_to_block_inner,
                            gf_struct_t const &gf_struct, double tau_split, double beta);

  frame_t make_u_frame(u_partial_t const &up);

  u_partial_t make_u_partial(frame_t const &u);

  frame_t make_frame(std::vector<long> const & shape);

  //
  frame_t make_frame(gf_struct_t const & gf_struct);

  u_partial_t operator*(u_partial_t const &l, u_partial_t const &r);

  u_partial_t apply_op_from_right(u_partial_t const &l, int lin_index, bool op_dag, atom_diag const & ad);

  // Calculate the Frobenius norm of the u_frame block diagonal matrix:
  double frobenius_norm(u_partial_t const &u_partial);
  double frobenius_norm(g_frame_t const &g_frame);

  // calculate the trace of the u_frame block diagonal matrix:
  double trace(frame_t const &u_frame);

  // calculate the relative distance of the two frames: ||l - r|| / max(||l||, ||r||)
  double relative_distance(frame_t const &l, frame_t const &r);

  // calculate the relative distance of two propagators
  // maximum relative distance of all frames
  double relative_distance(u_tau_t const &l, u_tau_t const &r);

  inline double trace(u_partial_t const &u_partial) { return trace(make_u_frame(u_partial)); }

  void print(frame_t const &u_frame, double factor);

  frame_t make_g_frame_from_l_and_r(atom_diag const &ad_imp, std::map<int, std::pair<int, int>> const &map_lin_idx_to_block_inner,
                                    gf_struct_t const &gf_struct, u_partial_t const &l, u_partial_t const &r);

} // namespace inchworm
