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
#include "./types.hpp"
#include "./params.hpp"
#include "./diagram/diagram.hpp"

#include <numeric>
#include <bitset>

//#include <inchworm/solver_core.hpp>
#include <triqs/gfs.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/utility/macros.hpp>
#include <triqs/arrays.hpp>

namespace inchworm {

  /// Necessary to use atom_diag block structure for the propagator.
  /** 
   * @param ad atom_diag of the system considered here.
   * @return The gf_struct necessary to call u_tau_t.
   */
  triqs::hilbert_space::gf_struct_t find_propagator_struct(atom_diag const &ad);

  struct propagator_frame {
    std::vector<matrix<dcomplex>> matrices; // The different matrices of the blocks.
    int acc_number;                         // Number of sample accumlated here

    // Constructor
    propagator_frame(triqs::atom_diag::atom_diag<false> const &ad);

    // Function to add them, and accumulate.
    // propagator_frame &operator+=(propagator_frame U_frame);

    // Function to add them, and accumulate.
    propagator_frame &operator+=(propagator_frame U_frame) {
      for (int bl = 0; bl < matrices.size(); bl++) matrices[bl] += U_frame.matrices[bl];
      acc_number++;
      return *this;
    }

    void assign(int bl, matrix<dcomplex> mat);

    // Set the values to zero
    void reset();

    // Calculate the Frobenius norm of the matrix
    double frobenius_norm();

    // Printing function.
    friend std::ostream &operator<<(std::ostream &out, propagator_frame const &U_frame);
  };

  void assign_frame_to_propagator(u_tau_t &U, propagator_frame const &U_frame, int frame);
  void assign_identity_to_propagator(u_tau_t &U, int frame);

  /// Function that calculate the product: U_frame = U(tau_0) op U(tau_1-tau_0) op U(tau_2-tau_1) op U(tau_3-tau_2) ... op U(tau-tau_n)
  /// where op is either c_dag or c operator, depending on the configuration
  /** 
   * @param U Full propagator calculated up until this point.
   * @param ad atom_diag of the system considered here.
   * @param diagram Configuration of the n operators (op) of the present Monte Carlo step.
   * @param tau Time of the propagator_frame calculated here. tau must be greater than any times
   * @param use_bare_U If true, calculate the same product using only the bare propagators. 
   * @return propagator_frame, at time tau, resulting from this product.
   */
  propagator_frame propagator_product(u_tau_t const &U, triqs::atom_diag::atom_diag<false> const &ad, time_diagram_t const &diagram, double tau,
                                      bool use_bare_U = false);

  /// If the user do not provide the propagator, use bare propagator instead.
  propagator_frame propagator_product(triqs::atom_diag::atom_diag<false> const &ad, time_diagram_t const &diagram, double tau);

} // namespace inchworm
