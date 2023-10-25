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

#include <triqs/test_tools/gfs.hpp>
#include <inchworm/diagram/inclusion_exclusion.hpp>
#include <inchworm/diagram/proper_enum.hpp>
#include <inchworm/diagram/hyb_matrix.hpp>
#include <inchworm/types.hpp>
#include <inchworm/diagram/print.hpp>

using namespace inchworm::diagram;
using namespace inchworm;

void compare_both_methods(std::vector<double> &tau1, std::vector<double> &tau2, std::vector<double> &split_times) {

  std::sort(tau1.begin(), tau1.end());
  std::sort(tau2.begin(), tau2.end());

  std::vector<fop_t> c, cdag;
  for (auto t : tau1) c.push_back({t, false, 0, 0, 0});
  for (auto t : tau2) cdag.push_back({t, true, 0, 0, 0});
  time_diagram_t diagram(c, cdag, split_times);
  auto hyb_mat = hyb_matrix_t(diagram);

  hyb_scalar_t value_inclus = inclusion_exclusion(diagram, hyb_mat);
}

std::vector<double> generate_random_vector(double beta, int n_tau) {
  //double beta = 50.0;
  std::vector<double> tau(n_tau);
  std::generate(tau.begin(), tau.end(), [beta]() mutable { return beta * (double)rand() / RAND_MAX; });
  return tau;
}

void print_vector(std::vector<double> const &v) {
  for (auto l : v) { std::printf("% 4.6f ", l); }
  std::printf("\n");
}

std::vector<double> generate_worst_case_scenario_vector(int n_tau, double first_value, double increment) {
  std::vector<double> tau(n_tau);
  std::generate(tau.begin(), tau.end(), [n = 0, first_value, increment]() mutable { return (first_value + (n++) * increment); });
  print_vector(tau);
  return tau;
}

int main(void) {
  int N         = 500000; //500000;
  int sp_max    = 1;
  int order_min = 5;  // order 0 is a special case that fails for now.
  int order_max = 10; // order 9 and above are quite slow
  double beta   = 1.0;

  for (int order = order_min; order <= order_max; order++) {
    std::printf("order = % d\n", order);
    for (int sp_number = 1; sp_number <= sp_max; sp_number++)
      for (int n = 0; n < N; n++) {

        std::vector<double> tau1        = generate_random_vector(beta, order);
        std::vector<double> tau2        = generate_random_vector(beta, order);
        std::vector<double> split_times = generate_random_vector(beta, sp_number);

        compare_both_methods(tau1, tau2, split_times);
      }
  }
  return 0;
}

/*
int main(void) {
  int N         = 1;//500000;
  int sp_max    = 1;
  int order_min = 5; 
  int order_max = 5; 
  double beta   = 1.0;

  for (int order = order_min; order <= order_max; order++)
  {
    std::printf("order = % d\n", order);
    for (int sp_number = 1; sp_number <= sp_max; sp_number++)
      for (int n = 0; n < N; n++) {

        std::vector<double> tau1        = generate_worst_case_scenario_vector(order, 0.0001, 0.0002);
        std::vector<double> tau2        = generate_worst_case_scenario_vector(order, 0.0000, 0.0002);
        std::vector<double> split_times = {0.000001};

        compare_both_methods(tau1, tau2, split_times);
      }
  }
  return 0;
}
*/
