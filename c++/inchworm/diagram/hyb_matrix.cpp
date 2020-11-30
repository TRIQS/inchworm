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
#include "./hyb_matrix.hpp"
#include "./utilities.hpp"

#include <iomanip>

namespace inchworm::diagram {

  using itertools::enumerate;

  hyb_matrix_t::hyb_matrix_t(time_diagram_t const &diagram) : diagram{diagram}, size{diagram.perturbation_order()}, mat{size, size} {

    auto hyb_function_dummy = [](double dtau) { return (1.0 / (0.8 * (dtau - 0.5))); };
    for (auto [i, d] : enumerate(diagram.d_list)) {
      for (auto [j, d_dag] : enumerate(diagram.d_dag_list)) {
        mat(i, j) = (0.5 + j - i) * hyb_function_dummy(d_dag.tau - d.tau);
        if (i < j) mat(i, j) = -mat(i, j);
      }
    }
  }

  hyb_matrix_t::hyb_matrix_t(time_diagram_t const &diagram, hyb_adaptor_t const &hyb_tau)
     : diagram{diagram}, size{diagram.perturbation_order()}, mat{size, size} {

    for (auto [i, d] : enumerate(diagram.d_list)) {
      for (auto [j, d_dag] : enumerate(diagram.d_dag_list)) { mat(i, j) = hyb_tau(d_dag, d); }
    }
  }

  void hyb_matrix_t::optimize_inclusion_exclusion() {
    if (smallest_segment == 4) {
      for (int k = 0; k < diagram.size() - 1; k++) {

        if (std::any_of(begin(diagram.split_points), end(diagram.split_points), [k](int i) { return i == k + 1; })) continue;

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
  }

  scalar_t hyb_matrix_t::det() const { return determinant(mat); }

  scalar_t hyb_matrix_t::extract_det(std::vector<int> const &list_of_indices) const {
    EXPECTS(list_of_indices.size() < 400); // I think that order 200 is safe enough.
    EXPECTS(list_of_indices.size() % 2 == 0 ); 
    
    int N = list_of_indices.size() / 2;
    //std::vector<int> list_of_d(100), list_of_d_dag(100);
    std::array<int, 400> list_of_d;
    std::array<int, 400> list_of_d_dag;
    
    int i = 0, j = 0;
    //std::printf("%d %d\n", i, j);
    for (int n=0; n<list_of_indices.size(); n++) {//auto idx : list_of_indices) {
      int idx = list_of_indices[n];
      if (diagram.op_list[idx].dag)
        list_of_d_dag[i++] = diagram.op_list[idx].order_index;
      else
        list_of_d[j++] = diagram.op_list[idx].order_index;
      //std::printf("%d %d\n", i, j);
    }

    //nda::matrix<scalar_t, 2, nda::C_layout, 'M', nda::sso<1000>> m(N, N);
    nda::matrix<scalar_t, nda::C_layout, nda::sso<1000>> m(N, N); // we might need to play with this number. For the test I did (order 5 to 10), it had minor speedup (~1%).
    //matrix_t m(N, N);

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
