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
#include "./segment.hpp"

#include <iomanip>

inline double one_fermion(double tau, double eps, double beta) {
  if (eps >= 0) {
    return -std::exp(-tau * eps) / (1. + std::exp(-beta * eps));
  } else {
    return -std::exp((beta - tau) * eps) / (1. + std::exp(beta * eps));
  }
}

namespace inchworm::diagram {

  using itertools::enumerate;

  hyb_matrix_t::hyb_matrix_t(time_diagram_t const &diagram) : diagram{diagram}, size{diagram.perturbation_order()}, mat{size, size} {

    // auto hyb_function_dummy = [](double dtau) { return (1.0 / (0.8 * (dtau - 0.5))); };
    // for (auto [i, d] : enumerate(diagram.d_list)) {
    //   for (auto [j, d_dag] : enumerate(diagram.d_dag_list)) {
    //     mat(i, j) = (0.5 + j - i) * hyb_function_dummy(d_dag.tau - d.tau);
    //     if (i < j) mat(i, j) = -mat(i, j);
    //   }
    // }
  }

  hyb_matrix_t::hyb_matrix_t(time_diagram_t const &diagram, hyb_tau_t const &Delta)
     : diagram{diagram}, size{diagram.perturbation_order()}, mat{size, size} {

    auto eval_Delta = [&Delta](fop_t const &cdag, fop_t const &c) -> hyb_scalar_t {
      if (cdag.bl != c.bl) return 0.;

      double dtau = cdag.tau - c.tau;
      if (dtau >= 0.) {
        return Delta[c.bl](dtau)(cdag.idx, c.idx);
      } else {
        return -Delta[c.bl](Delta[c.bl].mesh().beta() + dtau)(cdag.idx, c.idx);
      }

    };

    for (auto [i, d] : enumerate(diagram.d_list)) {
      for (auto [j, d_dag] : enumerate(diagram.d_dag_list)) { mat(i, j) = eval_Delta(d_dag, d); }
    }
  }

  hyb_matrix_t::hyb_matrix_t(time_diagram_t const &diagram, mat_t const &theta, vec_t const &eps, double beta, int n_bath)
     : diagram{diagram}, size{diagram.perturbation_order()}, mat{size, size} {

    auto eval_Delta = [&theta, &eps, &n_bath, &beta](fop_t const &cdag, fop_t const &c) -> hyb_scalar_t {
      if (cdag.bl != c.bl) return 0.;

      double dtau = cdag.tau - c.tau;
      int i = cdag.idx;
      int j = c.idx;
      double result=0.0;
      if (dtau >=0.){
        for (auto n : range(n_bath)) {
          result += theta(i, n) * dagger(theta)(n, j) * one_fermion(dtau, eps(n), beta);
        }
      }
      else{
       for (auto n : range(n_bath)) {
          result += -1.0 * theta(i, n) * dagger(theta)(n, j) * one_fermion(dtau + beta, eps(n), beta);
        }
      }
      return result;

    };

    for (auto [i, d] : enumerate(diagram.d_list)) {
      for (auto [j, d_dag] : enumerate(diagram.d_dag_list)) { mat(i, j) = eval_Delta(d_dag, d); }
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

  hyb_scalar_t hyb_matrix_t::det() const { return determinant(mat); }

  void hyb_matrix_t::print() const {
    std::printf("\nhybridization mat: \n");
    for (int i = 0; i < diagram.perturbation_order(); i++) {
      for (int j = 0; j < diagram.perturbation_order(); j++) { std::printf("% 2.7e ", mat(i, j)); }
      std::printf("\n");
    }
  }

  
} // namespace inchworm::diagram
