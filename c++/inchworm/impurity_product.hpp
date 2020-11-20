#pragma once
#include "./u_frame.hpp"
#include "./util.hpp"
#include "./types.hpp"
#include "./params.hpp"
#include "./diagram/diagram.hpp"

#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/utility/macros.hpp>
#include <triqs/arrays.hpp>

namespace inchworm {
  using time_diagram_t = diagram::time_diagram_t;

  /** 
   * Calculate the product: u_frame = u(tau_0) op u(tau_1-tau_0) op u(tau_2-tau_1) op u(tau_3-tau_2) ... op u(tau-tau_n)
   * where op is either c_dag or c operator, depending on the configuration
   *
   * @param ad atom_diag of the system considered here.
   * @param diagram configuration of the n operators (op) of the present monte carlo step.
   * @param tau_min The smallest time of the segment
   * @param tau_max The largest time of the segment
   * @param u_tau_p pointer to the full propagator (u_tau) calculated up until this point (0 < tau < tau_split). 
   * @return frame_t, at time tau_max, resulting from this product.
   */
  u_partial_t impurity_product(atom_diag const &ad, time_diagram_t const &diagram, double tau_max, double tau_min,
                               u_tau_t const *const u_tau_p = nullptr);

} // namespace inchworm
