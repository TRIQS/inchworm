#pragma once

#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/utility/macros.hpp>

#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/utility/time_pt.hpp>

#include <mpi/mpi.hpp>

#include <iostream>
#include <string>
#include <utility>

namespace inchworm {

  using namespace std::complex_literals; // Complex Unity 1i
  using namespace triqs::gfs;
  using namespace triqs::mesh;
  using namespace triqs::arrays;
  using namespace triqs::operators;
  using namespace triqs::hilbert_space;
  using namespace triqs::utility;
  using namespace h5;
  using namespace itertools;

  // Defined by Maxime, need comments
  // FIXME Add LOCAL_HAMILTONIAN_IS_COMPLEX ?
  // Carefully check type propagation
#ifdef HYBRIDISATION_IS_COMPLEX
  using scalar_t                            = dcomplex;
  static constexpr bool is_h_scalar_complex = true;
  using hyb_target_t                        = matrix_valued;
#else
  using scalar_t                            = double;
  static constexpr bool is_h_scalar_complex = false;
  using hyb_target_t                        = matrix_real_valued;
#endif
  using matrix_t = matrix<scalar_t>;

  /// The structure of the gf : block_idx -> pair of block_name and index list (int/string)
  using triqs::hilbert_space::gf_struct_t;

  /// Container type of the propagator
  using u_tau_t   = block_gf<imtime, hyb_target_t>;

  /// Container type of one-particle Green and Vertex functions in imaginary times
  using h_tau_t = block_gf<imtime, hyb_target_t>;

  /// Container type of one-particle Green and Vertex functions in imaginary times
  using g_tau_t = block_gf<imtime, hyb_target_t>;

  /// A view to a g_tau_t
  using g_tau_vt = g_tau_t::view_type;

  /// A const_view to a g_tau_t
  using g_tau_cvt = g_tau_t::const_view_type;

  /// Container type of one-particle Green and Vertex functions in Matsubara frequencies
  using g_iw_t = block_gf<imfreq, matrix_valued>;

  /// A view to a g_iw_t
  using g_iw_vt = g_iw_t::view_type;

  /// A const_view to a g_iw_t
  using g_iw_cvt = g_iw_t::const_view_type;

  /// Type of the Monte-Carlo weight. Either double or dcomplex
  //using scalar_t = scalar_t

  using atom_diag = triqs::atom_diag::atom_diag<is_h_scalar_complex>;

  using triqs::hilbert_space::gf_struct_t;
  //using triqs::utility::time_pt;
  //using op_t         = std::pair<time_pt, int>;
  using indices_type = triqs::operators::indices_t;

  // Declare some placeholders for the rest of the code. Use anonymous namespace for proper linkage
  // in this code, all variables with trailing _ are placeholders by convention.
  //constexpr triqs::clef::placeholder<0> iw_;
  //constexpr triqs::clef::placeholder<1> tau_;

} // namespace inchworm
