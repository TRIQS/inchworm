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

#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <numeric>

#include "diagram.hpp"
#include "../types.hpp"
#include "../mc/qmc_config_data.hpp"

namespace inchworm::diagram {
  //using hybridization_function_t = triqs::gfs::gf<triqs::gfs::imtime,triqs::gfs::matrix_real_valued>;
  using triqs::utility::enumerate;

  //scalar_t hyb_function(scalar_t dtau);

  struct hyb_matrix_t {

    using matrix_t = triqs::arrays::matrix<scalar_t>;

    matrix_t mat;
    time_diagram_t const &diagram;

    /// Constructors
    hyb_matrix_t(time_diagram_t const &diagram, std::function<scalar_t(double)>);
    hyb_matrix_t(time_diagram_t const &diagram, hyb_adaptor_t const &hyb_tau);

    /// Set the value of adjacent vertex to zero in the matrix. (segment of length 2 optimization)
    void optimize_inclusion_exclusion();

    /// Return determinant of full matrix
    scalar_t det();

    /// Extract a determinant of sub indices of the matrix
    scalar_t extract_det(std::vector<int> const &list_of_indices) const;

    void print();
  };
} // namespace inchworm::diagram
