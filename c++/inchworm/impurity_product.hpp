#pragma once

#include "types.hpp"
#include "diagram/diagram.hpp"

namespace inchworm {

  /**
   * Quickly check if the impurity trace defined by the diagram vanishes
   *
   * @param ad atom_diag object of the system under consideration
   * @param diagram The diagram configuration
   * @return True if it is zero, else false
   */
  bool has_zero_trace(atom_diag const &ad, diagram::time_diagram_t const & diagram);

  /** 
   * Calculate the operator product
   *
   *  u_frame = u(tau_max - tau_{n-1}) op_{n-1} u(tau_{n-1}-tau_{n-2}) op_{n-2} ... op_1 u(tau_1 - tau_0) op_0 u(tau_0 - tau_min)
   *
   * where (op_{n-1}, .., op_0) is a time-ordered list of c and c_dag operators defined by the configuration and (tau_max, tau_min),
   * such that tau_max < tau_{n-1} < .. < tau_0 < tau_min
   *
   * @param ad atom_diag object of the system under consideration
   * @param diagram The diagram configuration
   * @param tau_min The smallest time of the segment
   * @param tau_max The largest time of the segment
   * @param u_tau_p Pointer to the full propagator. Must be initialized for all 0 < tau < tau_max - tau_min
   * @return The operator product
   */
  u_partial_t impurity_product(atom_diag const &ad, diagram::time_diagram_t const &diagram, double tau_max, double tau_min,
                               u_tau_t const *const u_tau_p = nullptr);

} // namespace inchworm
