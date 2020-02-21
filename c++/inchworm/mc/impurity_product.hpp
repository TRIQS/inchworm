/*******************************************************************************
 *
 * inchworm: A TRIQS based impurity solver
 *
 * Copyright (c) 2019 The Simons foundation
 *   authors: Nils Wentzell
 *
 * inchworm is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * inchworm is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * inchworm. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once
#include "../types.hpp"
#include "../params.hpp"
#include "../diagram/diagram.hpp"

#include <numeric>
#include <bitset>

//#include <inchworm/solver_core.hpp>
#include <triqs/gfs.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/utility/macros.hpp>
#include <triqs/arrays.hpp>

namespace inchworm {
  using time_diagram_t = diagram::time_diagram_t;

  /// Necessary to use atom_diag block structure for the propagator.
  /** 
   * @param ad atom_diag of the system considered here.
   * @return The gf_struct necessary to call u_tau_t.
   */
  gf_struct_t find_propagator_struct(atom_diag const &ad);

  u_frame_t make_zero_propagator_frame(atom_diag const &ad);
  u_frame_t make_bare_propagator_frame(atom_diag const &ad, double tau);

  double frobenius_norm(u_frame_t const &u_frame);
  double trace(u_frame_t const &u_frame);

  //void init_propagator_frame(u_tau_t &U, u_frame_t const &u_frame, int frame); // faire un constructeur (struct)

  u_tau_t make_propagator(atom_diag const &h_diag, int n_tau);

  /// Function that calculate the product: u_frame = U(tau_0) op U(tau_1-tau_0) op U(tau_2-tau_1) op U(tau_3-tau_2) ... op U(tau-tau_n)
  /// where op is either c_dag or c operator, depending on the configuration
  /** 
   * @param ad atom_diag of the system considered here.
   * @param diagram Configuration of the n operators (op) of the present Monte Carlo step.
   * @param tau Time of the u_frame_t calculated here. tau must be greater than any times
   * @param u_tau Full propagator calculated up until this point.
   * @return u_frame_t, at time tau, resulting from this product.
   */
  u_frame_t propagator_product(atom_diag const &ad, time_diagram_t const &diagram, double tau, u_tau_t const *const u_tau_p = nullptr);

  inline std::ostream &operator<<(std::ostream &out, u_frame_t const &u_frame) {
    out << "propagator_frame (size: " << u_frame.size() << ")\n";
    for (int bl = 0; bl < u_frame.size(); bl++) { out << u_frame[bl] << "\n"; }
    return out;
  }

  //constexpr int MAX_ORDER = 10;
  struct single_step_results_t {
    double average_sign = 0.0;
    double average_k    = 0.0;
    u_frame_t u_frame;
    u_frame_t zero_frame;
    //std::vector<double> u_expansion_order;
    single_step_results_t(atom_diag const &h_diag) { //: u_expansion_order(MAX_ORDER, 0) {
      u_frame = make_zero_propagator_frame(h_diag);
      zero_frame = u_frame;
    };
  };

} // namespace inchworm
