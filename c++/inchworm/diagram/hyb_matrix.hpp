/*******************************************************************************
 *
 * inchworm: A TRIQS based impurity solver
 *
 * Copyright (c) 2019 The Simons foundation
 *   authors: Maxime Charlebois
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

#include "diagram.hpp"

#include "../types.hpp"
#include "../config.hpp"

namespace inchworm::diagram {

  struct hyb_matrix_t {

    /// Construct hybridization matrix using dummy function (for test purposes only)
    hyb_matrix_t(time_diagram_t const &diagram);

    /// Construct hybridization matrix using a hyb_adaptor_t
    hyb_matrix_t(time_diagram_t const &diagram, hyb_tau_t const &Delta);

    /// Optimization for segments of length 2:
    ///   Set the value of adjacent vertices to zero in the matrix.
    void optimize_inclusion_exclusion();

    /// Return determinant of full matrix
    hyb_scalar_t det() const;

    /// Calculate the determinant of the submatrix defined by list_of_indices
    template <typename R> hyb_scalar_t extract_det(R const &list_of_indices) const;

    void print() const;

    /// The associated diagram
    time_diagram_t const &diagram;

    /// The linear size of the matrix
    int size;

    /// The matrix
    matrix<hyb_scalar_t> mat;
  };

  template <typename R> hyb_scalar_t hyb_matrix_t::extract_det(R const &list_of_indices) const {
    EXPECTS(list_of_indices.size() % 2 == 0);

    // Creation of C-Style arrays slightly more performant than sso nda::array
    int N = list_of_indices.size() / 2;
    int list_of_d[N], list_of_d_dag[N];

    for (int i = 0, j = 0; auto idx : list_of_indices) {
      if (diagram.op_list[idx].dag)
        list_of_d_dag[i++] = diagram.op_list[idx].order_index;
      else
        list_of_d[j++] = diagram.op_list[idx].order_index;
    }

    nda::matrix<hyb_scalar_t, nda::C_layout, nda::sso<1000>> m(N, N);
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) { m(i, j) = mat(list_of_d[i], list_of_d_dag[j]); }
    }

    return nda::linalg::det_in_place(m);
  }

} // namespace inchworm::diagram
