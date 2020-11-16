#pragma once
#include "./u_frame.hpp"
#include "./util.hpp"
#include "./types.hpp"
#include "./params.hpp"
#include "./diagram/diagram.hpp"

//#include <inchworm/solver_core.hpp>
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/utility/macros.hpp>
#include <triqs/arrays.hpp>

namespace inchworm {
  using time_diagram_t = diagram::time_diagram_t;

  // FIXME: these three functions should be put in a seperated file
  // TODO: Move this file out of the mc directory

  // Make an empty propagator (green function) with the same structure as the one in atom_diag:
  u_tau_t make_propagator(atom_diag const &ad_imp, double beta, int n_tau);

  // Make an exact diagonalization propagator U = Trace_B [exp(-H_bath *(beta-tau)) exp(-H_tot*tau)  ]  /  Trace_B [ exp(-H_bath*beta) ]
  u_tau_t make_ED_propagator(atom_diag const &ad_tot, atom_diag const &ad_imp, atom_diag const &ad_bath, double beta, int n_tau);

  // zeroth order is :
  // U_0(tau) = exp(-H_imp*tau) for cthyb
  // U_0(tau) = U(tau-tau_split) U(tau_split) for inchworm
  frame_t make_zeroth_order_frame(atom_diag const &ad, double tau, double tau_split = 0.0, u_tau_t const *const u_tau_p = nullptr);

  /// function that calculate the product: u_frame = u(tau_0) op u(tau_1-tau_0) op u(tau_2-tau_1) op u(tau_3-tau_2) ... op u(tau-tau_n)
  /// where op is either c_dag or c operator, depending on the configuration
  /** 
   * @param ad atom_diag of the system considered here.
   * @param diagram configuration of the n operators (op) of the present monte carlo step.
   * @param tau_min The smallest time of the segment
   * @param tau_max The largest time of the segment
   * @param u_tau_p pointer to the full propagator (u_tau) calculated up until this point (0 < tau < tau_split). 
   * @return frame_t, at time tau_max, resulting from this product.
   */
  u_partial_t impurity_product(atom_diag const &ad, time_diagram_t const &diagram, double tau_min, double tau_max,
                               u_tau_t const *const u_tau_p = nullptr);

} // namespace inchworm
