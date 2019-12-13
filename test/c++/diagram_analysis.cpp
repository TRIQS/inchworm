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
#include <inchworm/diagram/proper_enum.hpp>
#include <inchworm/diagram/inclusion_exclusion.hpp>

void compare_both_methods(std::vector<double> &tau1, std::vector<double> &tau2, std::vector<double> &split_times) {

  std::sort(tau1.begin(), tau1.end());
  std::sort(tau2.begin(), tau2.end());

  std::vector<time_and_orbital_t> c, cdag;
  for (auto t : tau1) c.push_back({t, 0});
  for (auto t : tau2) cdag.push_back({t, 0});
  time_diagram_t diagram(c, cdag, split_times);

  //int N_proper = find_proper_diagrams(diagram);
  hybridization_scalar_t value_proper = proper_enum(diagram);
  if constexpr (verbose) std::printf("c_k = % 4.6f\n", value_proper);

  hybridization_scalar_t value_inclus = inclusion_exclusion(diagram);
  if constexpr (verbose) std::printf("c_k = % 4.6f\n\n\n", value_inclus);

  if constexpr (verbose) std::printf("proper-enum         c_k = % 4.6f\n", value_proper);
  if constexpr (verbose) std::printf("inclusion-exclusion c_k = % 4.6f\n\n", value_inclus);

  EXPECT_NEAR(value_proper, value_inclus,
              std::max(1e-8, std::abs(1e-8 * value_proper))); //note: according to my random tests, 1e-9 was too strick in some extreme cases
}

std::vector<double> generate_random_vector(double beta, int n_tau) {
  //double beta = 50.0;
  std::vector<double> tau(n_tau);
  std::generate(tau.begin(), tau.end(), [beta]() mutable { return beta * (double)rand() / RAND_MAX; });
  return tau;
}

TEST(inchworm, benchmark_both_random) {
  int N         = 20;
  int sp_max    = 4;
  int order_min = 2; //9;   // order 0 and 1 are special case that fails for now.
  int order_max = 8; //10;
  double beta   = 1.0;

  for (int order = order_min; order <= order_max; order++)
    for (int sp_number = 1; sp_number < sp_max; sp_number++)
      for (int n = 0; n < N; n++) {

        std::vector<double> tau1        = generate_random_vector(beta, order);
        std::vector<double> tau2        = generate_random_vector(beta, order);
        std::vector<double> split_times = generate_random_vector(beta, sp_number);

        compare_both_methods(tau1, tau2, split_times);
      }
}

/*
TEST(inchworm, benchmark_both1) {

  std::vector<double> tau1 = generate_random_vector(1.0, 7);
  std::vector<double> tau2 = generate_random_vector(1.0, 7);
  std::vector<double> split_times = generate_random_vector(1.0, 1);

  compare_both_methods(tau1, tau2, split_times);
}

TEST(inchworm, benchmark_both1) {
  std::vector<double> tau1        = {0.1, 0.3, 0.4344, 0.5, 0.75, 0.9};
  std::vector<double> tau2        = {0, 0.21, 0.3452, 0.45, 0.69, 0.81};
  std::vector<double> split_times = {0.8, 0.55, 0.34};
  compare_both_methods(tau1, tau2, split_times);
}

TEST(inchworm, benchmark_both2) {
  std::vector<double> tau1        = {1.0, 2.0, 3.0, 4.5};
  std::vector<double> tau2        = {1.2, 1.8, 4.6, 4.7};
  std::vector<double> split_times = {1.9};
  compare_both_methods(tau1, tau2, split_times);
}

//long test, comment if verbose > 1
TEST(inchworm, benchmark_both3) {
  std::vector<double> tau1        = {0.1, 0.5, 0.6, 0.89, 1.3, 1.45, 1.78, 4.2};
  std::vector<double> tau2        = {0.0, 0.11, 0.51, 1.2, 1.8, 2.2, 2.3, 4.6};
  std::vector<double> split_times = {0.99};
  compare_both_methods(tau1, tau2, split_times);
  split_times.push_back(0.2);
  compare_both_methods(tau1, tau2, split_times);
  split_times.push_back(3.7);
  compare_both_methods(tau1, tau2, split_times);
}
//

//long test, comment if verbose > 1
TEST(inchworm, benchmark_both4) {
  std::vector<double> tau1 = {0.1, 0.2, 0.3, 0.5, 1.3, 1.4, 1.5, 3.6};
  std::vector<double> tau2 = {0.0, 0.53, 0.54, 1.2, 1.8, 2.2, 2.3, 4.6};
  std::vector<double> split_times = {0.99};
  compare_both_methods(tau1, tau2, split_times);
  split_times.push_back(0.21);
  compare_both_methods(tau1, tau2, split_times);
  split_times.push_back(3.7);
  compare_both_methods(tau1, tau2, split_times);
}


TEST(inchworm, inclusion_exclusion_big_order1) {
  //std::vector<double> tau1={0.4344,0.1,0.3,0.5,0.75,0.9,0.55,0.566,0.33,.4959594,.4494929,.12349512,.62343,0.123412444,0.2134444,.99949941,1.1,1.23,1.45,2.3,1.6,4.3,1.222,2.98,3.1244};
  //std::vector<double> tau2={0,0.3452,0.21,0.45,0.69,0.81,0.998,0.122,0.833,0.4934,.210342134,.210343,.02134,.0030404,.02142430,0.1111,1.2,1.3,1.4,3.4,2.34,3.11,2.9,1.99,3.098};

  std::vector<double> tau1 = {0.4344, 0.1,      0.3,      0.5,       0.75,   0.9,         0.55,      0.566,
                              0.33,   .4959594, .4494929, .12349512, .62343, 0.123412444, 0.2134444, .99949941};
  std::vector<double> tau2 = {0,     0.3452, 0.21,       0.45,    0.69,   0.81,     0.998,     0.122,
                              0.833, 0.4934, .210342134, .210343, .02134, .0030404, .02142430, 0.1111};
  std::vector<double> split_times = {0.82};

  std::sort(tau1.begin(), tau1.end());
  std::sort(tau2.begin(), tau2.end());

  std::vector<time_and_orbital_t> c, cdag;
  for (auto t : tau1) c.push_back({t, 0});
  for (auto t : tau2) cdag.push_back({t, 0});
  time_diagram_t diagram(c, cdag, split_times);

  hybridization_scalar_t value = inclusion_exclusion(diagram);
  std::printf("c_k = % 4.6f\n\n", value);

  EXPECT_NEAR(-8458508.82790539, value, 1e-7);
}
//*/
/*
TEST(inchworm, inclusion_exclusion_huge_order1) {
  std::vector<double> tau1={0.4344,0.1,0.3,0.5,0.92882,0.75,0.9,0.55,0.566,0.33,.4959594,.4494929,.12349512,.62343,0.123412444,0.2134444,.99949941,1.1,1.23,1.45,2.3,1.6,4.3,1.222,2.98,3.1244};
  std::vector<double> tau2={0,0.3452,0.21,0.45,0.9329,0.69,0.81,0.998,0.122,0.833,0.4934,.210342134,.210343,.02134,.0030404,.02142430,0.1111,1.2,1.3,1.4,3.4,2.34,3.11,2.9,1.99,3.098};

  //std::vector<double> tau1={0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1};
  //std::vector<double> tau2={0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0,5.2};
  //std::vector<double> tau1={0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1};
  //std::vector<double> tau2={0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2};

  //std::vector<double> tau1 = {0.4344, 0.1,      0.3,      0.5,       0.75,   0.9,         0.55,      0.566,
  //                            0.33,   .4959594, .4494929, .12349512, .62343, 0.123412444, 0.2134444, .99949941};
  //std::vector<double> tau2 = {0,     0.3452, 0.21,       0.45,    0.69,   0.81,     0.998,     0.122,
  //                            0.833, 0.4934, .210342134, .210343, .02134, .0030404, .02142430, 0.1111};
  std::vector<double> split_times = {0.82};

  std::sort(tau1.begin(), tau1.end());
  std::sort(tau2.begin(), tau2.end());

  std::vector<time_and_orbital_t> c, cdag;
  for (auto t : tau1) c.push_back({t, 0});
  for (auto t : tau2) cdag.push_back({t, 0});
  time_diagram_t diagram(c, cdag, split_times);

  hybridization_scalar_t value = inclusion_exclusion(diagram);
  std::printf("c_k = % 4.6f\n\n", value);

  //EXPECT_NEAR(-8458508.82790539, value, 1e-7);
}
//*/
MAKE_MAIN
