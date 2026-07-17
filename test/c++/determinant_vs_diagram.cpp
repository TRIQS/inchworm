// Copyright (c) 2020--present, The Simons Foundation
// This file is part of inchworm and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.


#include <triqs/test_tools/gfs.hpp>
#include <inchworm/diagram/inclusion_exclusion.hpp>
#include <inchworm/diagram/proper_enum.hpp>
#include <inchworm/diagram/hyb_matrix.hpp>
#include <inchworm/types.hpp>

using namespace inchworm::diagram;
using namespace inchworm;

void compare_proper_with_det(std::vector<double> &tau1, std::vector<double> &tau2) {

  std::sort(tau1.begin(), tau1.end());
  std::sort(tau2.begin(), tau2.end());

  std::vector<fop_t> c, cdag;
  for (auto t : tau1) c.push_back({t, false, 0, 0, 0});
  for (auto t : tau2) cdag.push_back({t, false, 0, 0, 0});
  time_diagram_t diagram(c, cdag, {});
  auto hyb_mat = hyb_matrix_t(diagram);

  hyb_scalar_t value_det = hyb_mat.det();
  hyb_mat.print();
  //int N_proper = find_proper_diagrams(diagram);
  hyb_scalar_t value_proper = full_enum(diagram, hyb_mat, true /*verbose*/);
  std::printf("full-enum c_k   = % 4.6e\n", value_proper);
  std::printf("determinant c_k = % 4.6e\n\n", value_det);

  EXPECT_COMPLEX_NEAR(value_proper, value_det,
                      std::max(1e-8, std::abs(1e-8 * value_proper))); //note: according to my random tests, 1e-9 was too strick in some extreme cases
}

std::vector<double> generate_random_vector(double beta, int n_tau) {
  //double beta = 50.0;
  std::vector<double> tau(n_tau);
  std::generate(tau.begin(), tau.end(), [beta]() mutable { return beta * (double)rand() / RAND_MAX; });
  return tau;
}

TEST(inchworm, benchmark_both2) {
  std::vector<double> tau1 = {0.1};
  std::vector<double> tau2 = {0.05};
  compare_proper_with_det(tau1, tau2);
}

TEST(inchworm, benchmark_both_random) {
  int N = 4;
  //int sp_max    = 4;
  int order_min = 2; // order 0 is a special case that fails for now.
  int order_max = 5; // order 9 and above are quite slow
  double beta   = 1.0;

  for (int order = order_min; order <= order_max; order++)
    for (int n = 0; n < N; n++) {

      std::vector<double> tau1 = generate_random_vector(beta, order);
      std::vector<double> tau2 = generate_random_vector(beta, order);

      compare_proper_with_det(tau1, tau2);
    }
}

MAKE_MAIN
