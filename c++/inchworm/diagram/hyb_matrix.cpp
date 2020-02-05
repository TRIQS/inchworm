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
  //scalar_t hyb_function_dummy(double dtau) { // to be changed in the future
  //  return (1.0 / (0.1 * dtau - 0.5));
  //}

  /*
  hyb_matrix_t::hyb_matrix_t(time_diagram_t const &diagram)
     : mat(diagram.perturbation_order(), diagram.perturbation_order()), diagram{diagram} {

    for (auto [i, c] : enumerate(diagram.c_list))
      for (auto [j, cdag] : enumerate(diagram.cdag_list)) {

        double dtau = cdag.tau - c.tau;
        mat(i, j)   = hyb_function(dtau);
      }
  }
*/

  hyb_matrix_t::hyb_matrix_t(time_diagram_t const &diagram, std::function<scalar_t(double)> hyb_function)
     //hyb_matrix_t::hyb_matrix_t(time_diagram_t const &diagram)
     : mat(diagram.perturbation_order(), diagram.perturbation_order()), diagram{diagram} {

    for (auto [i, c] : enumerate(diagram.c_list))
      for (auto [j, cdag] : enumerate(diagram.cdag_list)) {
        //for(auto const & cdag : diagram.cdag_list){

        double dtau = cdag.tau - c.tau;
        mat(i, j)   = hyb_function(dtau);
        //put this below in hyb_funciton at some point:
        //        if(dtau>=0) mat(i,j) = hyb( dtau )( cdag.orb, c.orb );
        //        else mat(i,j) = -hyb( hyb.mesh().domain().beta + dtau )( cdag.orb, c.orb );
      }
  }

  hyb_matrix_t::hyb_matrix_t(time_diagram_t const &diagram, hyb_adaptor_t const &hyb_tau)
     : mat(diagram.perturbation_order(), diagram.perturbation_order()), diagram{diagram} {

    for (auto [i, c] : enumerate(diagram.c_list))
      for (auto [j, cdag] : enumerate(diagram.cdag_list)) {

        //double dtau = cdag.tau - c.tau;
        mat(i, j) = hyb_tau(c.tau, c.linear_index, cdag.tau, cdag.linear_index);
      }
  }

  void hyb_matrix_t::optimize_inclusion_exclusion() {
    if (smallest_segment == 4) {
      for (int k = 0; k < diagram.op_list.size() - 1; k++) {

        if (std::any_of(begin(diagram.split_points), end(diagram.split_points), [k](int l) { return l == k + 1; })) continue;

        if (not diagram.op_list[k].dag and diagram.op_list[k + 1].dag) { // segment of length 2
          int it          = diagram.op_list[k + 1].order_index;
          int it_dag      = diagram.op_list[k].order_index;
          mat(it, it_dag) = 0.;
          //if constexpr (verbose) std::printf("order_indices  %d %d\n", diagram.op_list[k + 1].order_index, diagram.op_list[k].order_index);
        } else if (diagram.op_list[k].dag and not diagram.op_list[k + 1].dag) { // segment on length 2
          int it          = diagram.op_list[k].order_index;
          int it_dag      = diagram.op_list[k + 1].order_index;
          mat(it, it_dag) = 0.;
          //if constexpr (verbose) std::printf("order_indices   %d %d\n", diagram.op_list[k].order_index, diagram.op_list[k+1].order_index);
        }
      }
    }
    if constexpr (verbose > 1) print();
  }

  scalar_t hyb_matrix_t::det() { return determinant(mat); }

  scalar_t hyb_matrix_t::extract_det(std::vector<int> const &list_of_indices) const {

    int N = list_of_indices.size() / 2;
    EXPECTS(list_of_indices.size() % 2 == 0);
    matrix_t m(N, N);

    std::vector<int> list_of_c(N), list_of_cdag(N); // to refactor

    int i = 0, j = 0;
    for (auto idx : list_of_indices) {
      if (diagram.op_list[idx].dag)
        list_of_cdag[i++] = diagram.op_list[idx].order_index;
      else
        list_of_c[j++] = diagram.op_list[idx].order_index;
    }

    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) { m(i, j) = mat(list_of_cdag[i], list_of_c[j]); }
    }

    return determinant(m);
  }

  void hyb_matrix_t::print() {
    std::printf("\nhybridization mat: \n");
    for (int i = 0; i < diagram.perturbation_order(); i++) {
      for (int j = 0; j < diagram.perturbation_order(); j++) { std::printf("% 2.5f ", mat(i, j)); }
      std::printf("\n");
    }
    std::cout << std::setprecision(10) << mat;
    std::printf("\ndet: %e\n", det());
  }
} // namespace inchworm::diagram
