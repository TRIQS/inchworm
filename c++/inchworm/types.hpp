#pragma once

#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/atom_diag/functions.hpp>
#include <triqs/stat/lin_binning.hpp>
#include <triqs/stat/log_binning.hpp>

#include <nda/nda.hpp>
#include <nda/macros.hpp>
#include <nda/clef/literals.hpp>

#include <iomanip>
#include <concepts>

namespace inchworm {

  using namespace std::complex_literals; // Complex Unity 1i
  using namespace triqs::gfs;
  using namespace triqs::mesh;
  using namespace nda;
  using namespace triqs::operators;
  using namespace triqs::hilbert_space;
  using namespace triqs::utility;
  using namespace triqs::stat;
  using namespace h5;
  using namespace itertools;
  using namespace nda::clef::literals;

  /// The value type of the hybridization function
#ifdef HYBRIDIZATION_IS_COMPLEX
  using hyb_scalar_t = dcomplex;
  using hyb_target_t = matrix_valued;
#else
  using hyb_scalar_t = double;
  using hyb_target_t = matrix_real_valued;
#endif

  /// The value type of the impurity Hamiltonian
#ifdef IMPURITY_HAMILTONIAN_IS_COMPLEX
  using h_scalar_t = dcomplex;
#else
  using h_scalar_t = double;
#endif

  /// The target type of the impurity Green function and propagator
#if defined(HYBRIDIZATION_IS_COMPLEX) || defined(IMPURITY_HAMILTONIAN_IS_COMPLEX)
  using target_t = matrix_valued;
#else
  using target_t = matrix_real_valued;
#endif

  // Combined scalar_t: complex if either hyb or himp are complex
  using scalar_t = decltype(h_scalar_t{} * hyb_scalar_t{});
  using matrix_t = nda::matrix<scalar_t>;

  enum class MODE { PROPAGATOR, GREENFUNCTION };

  /// Short-hand for SSO vector type
  template <typename value_t> using sso_vector = nda::vector<value_t, nda::sso<100>>;

  /// The frame of a Green function or Propagator
  using frame_t = nda::array<matrix_t, 1>;

  // FIXME??? using u_partial_t = nda::array<std::pair<int, matrix_t>, 1>;
  using u_partial_t = std::vector<std::pair<int, matrix_t>>;

  /// The structure of the gf : block_idx -> pair of block_name and index list (int/string)
  using triqs::hilbert_space::gf_struct_t;

  // The many body operator type
  using many_body_op_t = triqs::operators::many_body_operator_generic<scalar_t>;

  /// Container type of the propagator
  using u_tau_t = block_gf<imtime, target_t>;

  /// Container type of one-particle Green and Vertex functions in imaginary times
  using hyb_tau_t = block_gf<imtime, hyb_target_t>;

  /// Container type of one-particle Green and Vertex functions in imaginary times
  using g_tau_t = block_gf<imtime, target_t>;

  /// Container type of one-particle Green and Vertex functions in Matsubara frequencies
  using g_iw_t = block_gf<imfreq, target_t>;

  /// The atom diag type to use
  using atom_diag = triqs::atom_diag::atom_diag<std::is_same_v<scalar_t, dcomplex>>;
  using triqs::hilbert_space::gf_struct_t;
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

    /// The left and right width of the insertion-move tau distribution
    double left_width  = 100.0;
    double right_width = 100.0;

    bool operator==(fop_t const &) const = default;

    friend inline bool operator<(fop_t const &o1, fop_t const &o2) { return o1.tau < o2.tau; }
  };

  inline std::ostream &operator<<(std::ostream &os, fop_t const &op) {
    os << std::setprecision(4);
    os << "c";
    if (op.dag) os << "_dag";
    os << "[tau: " << op.tau << ", bl: " << op.bl << ", idx: " << op.idx << ", lw: " << op.left_width << ", rw: " << op.right_width << "]";
    return os;
  }

} // namespace inchworm
