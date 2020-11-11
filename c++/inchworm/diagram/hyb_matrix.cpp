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
#include "hyb_matrix.hpp"
#include <iomanip>

namespace inchworm::diagram {
  scalar_t hyb_function_dummy(double dtau) { // for tests purpose only
    return (1.0 / (0.8 * (dtau - 0.5)));
  }

  // For stand alone tests only:
  hyb_matrix_t::hyb_matrix_t(time_diagram_t const &diagram) : mat(diagram.perturbation_order(), diagram.perturbation_order()), diagram{diagram} {
    int N = 0;
    for (auto [i, d] : enumerate(diagram.d_list)) {
      N++;
      for (auto [j, d_dag] : enumerate(diagram.d_dag_list)) {

        double dtau = d_dag.tau - d.tau;
        mat(i, j)   = (0.5 + j - i) * hyb_function_dummy(dtau);
        if (i < j) mat(i, j) = -mat(i, j);
      }
    }
    size = N;
  }

  // Extract hyb from adaptor and build hyb matrix:
  hyb_matrix_t::hyb_matrix_t(time_diagram_t const &diagram, hyb_adaptor_t const &hyb_tau)
     : mat(diagram.perturbation_order(), diagram.perturbation_order()), diagram{diagram} {

    for (auto [i, d] : enumerate(diagram.d_list)) {
      for (auto [j, d_dag] : enumerate(diagram.d_dag_list)) { mat(i, j) = hyb_tau(d_dag, d); }
    }
    size = diagram.d_list.size();
  }

  // optimization, set to zero components of the matrix corresponding to segment of length 2.
  void hyb_matrix_t::optimize_inclusion_exclusion() {
    if (smallest_segment == 4) {
      for (int k = 0; k < diagram.op_list.size() - 1; k++) {

        if (std::any_of(begin(diagram.split_points), end(diagram.split_points), [k](int l) { return l == k + 1; })) continue;

        if (not diagram.op_list[k].dag and diagram.op_list[k + 1].dag) { // segment of length 2
          int it          = diagram.op_list[k + 1].order_index;
          int it_dag      = diagram.op_list[k].order_index;
          mat(it_dag, it) = 0.;
        } else if (diagram.op_list[k].dag and not diagram.op_list[k + 1].dag) { // segment on length 2
          int it          = diagram.op_list[k].order_index;
          int it_dag      = diagram.op_list[k + 1].order_index;
          mat(it_dag, it) = 0.;
        }
      }
    }
    //if constexpr (verbose > 1) print();
  }

  scalar_t hyb_matrix_t::det() {
    /*    if (size == 0)
      return 1.0;
    else if (size == 1)
      return mat(0, 0);
    else if (size == 2)
      return (mat(0,0)*mat(1,1)-mat(1,0)*mat(0,1));
    else */
    return determinant(mat);
  }

  // Extract the submatrix composed of indices "list_of_indices" and compute the determinant of this submatrix:
  scalar_t hyb_matrix_t::extract_det(std::vector<int> const &list_of_indices) const {

    int N = list_of_indices.size() / 2;
    EXPECTS(list_of_indices.size() % 2 == 0);
    matrix_t m(N, N);

    std::vector<int> list_of_d(N), list_of_d_dag(N); // to refactor

    int i = 0, j = 0;
    for (auto idx : list_of_indices) {
      if (diagram.op_list[idx].dag)
        list_of_d_dag[i++] = diagram.op_list[idx].order_index;
      else
        list_of_d[j++] = diagram.op_list[idx].order_index;
    }

    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) { m(i, j) = mat(list_of_d[i], list_of_d_dag[j]); }
    }

    return determinant(m);
  }

  void hyb_matrix_t::print() const {
    std::printf("\nhybridization mat: \n");
    for (int i = 0; i < diagram.perturbation_order(); i++) {
      for (int j = 0; j < diagram.perturbation_order(); j++) { std::printf("% 2.7e ", mat(i, j)); }
      std::printf("\n");
    }
  }
} // namespace inchworm::diagram
