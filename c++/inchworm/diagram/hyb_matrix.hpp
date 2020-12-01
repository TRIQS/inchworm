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
#include "../qmc_config_data.hpp"

#define ORDER_MAX 400

namespace inchworm::diagram {

  struct hyb_matrix_t {

    /// Construct hybridization matrix using dummy function (for test purposes only)
    hyb_matrix_t(time_diagram_t const &diagram);

    /// Construct hybridization matrix using a hyb_adaptor_t
    hyb_matrix_t(time_diagram_t const &diagram, hyb_adaptor_t const &hyb_tau);

    /// Optimization for segments of length 2:
    ///   Set the value of adjacent vertices to zero in the matrix.
    void optimize_inclusion_exclusion();

    /// Return determinant of full matrix
    scalar_t det() const;

    /// Calculate the determinant of the submatrix defined by list_of_indices
    scalar_t extract_det(std::vector<int> const &list_of_indices) const;

    void print() const;

    /// The associated diagram
    time_diagram_t const &diagram;

    /// The linear size of the matrix
    int size;

    /// The matrix
    matrix_t mat;
  };
} // namespace inchworm::diagram
