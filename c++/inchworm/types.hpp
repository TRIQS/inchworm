#pragma once

#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/utility/macros.hpp>

#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/utility/time_pt.hpp>

#include <nda/nda.hpp>

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

  /// The frame of a Green function or Propagator
  using frame_t     = nda::array<matrix_t, 1>;
  using u_partial_t = std::vector<std::pair<int, matrix_t>>;

  /// The structure of the gf : block_idx -> pair of block_name and index list (int/string)
  using triqs::hilbert_space::gf_struct_t;

  /// Container type of the propagator
  using u_tau_t = block_gf<imtime, hyb_target_t>;

  /// A view to a u_tau_t
  using u_tau_vt = u_tau_t::view_type;

  /// A const view to a u_tau_t
  using u_tau_cvt = u_tau_t::const_view_type;

  /// Container type of one-particle Green and Vertex functions in imaginary times
  using h_tau_t = block_gf<imtime, hyb_target_t>;

  /// A view to a h_tau_t
  using h_tau_vt = h_tau_t::view_type;

  /// A const view to a h_tau_t
  using h_tau_cvt = h_tau_t::const_view_type;

  /// Container type of one-particle Green and Vertex functions in imaginary times
  using g_tau_t = block_gf<imtime, hyb_target_t>;

  /// A view to a g_tau_t
  using g_tau_vt = g_tau_t::view_type;

  /// A const view to a g_tau_t
  using g_tau_cvt = g_tau_t::const_view_type;

  /// Container type of one-particle Green and Vertex functions in Matsubara frequencies
  using g_iw_t = block_gf<imfreq, matrix_valued>;

  /// A view to a g_iw_t
  using g_iw_vt = g_iw_t::view_type;

  /// A const_view to a g_iw_t
  using g_iw_cvt = g_iw_t::const_view_type;

  /// The atom diag type to use
  using atom_diag = triqs::atom_diag::atom_diag<is_h_scalar_complex>;

  using indices_type = triqs::operators::indices_t;

  /**
   * Type representing a creation or annihilation operator at a fixed time
   */
  struct fop_t {

    /// The imaginary time // FIXME Or time_pt ?
    double tau;

    /// C (false) or Cdag (true) // FIXME Can we make this compile-time?
    bool dag;

    /// The linear index in the fundamental_operator_set
    long linear_index; // FIXME Can we get rid of this?

    /// The block index
    long bl;

    /// The orbital (or non-block) index
    long idx;

    friend inline bool operator<(fop_t const &o1, fop_t const &o2) { return o1.tau < o2.tau; }
  };

} // namespace inchworm
