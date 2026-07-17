// Copyright (c) 2021--present, The Simons Foundation
// This file is part of inchworm and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.


#include <nda/gtest_tools.hpp>

#include <inchworm/torus.hpp>
#include <inchworm/config.hpp>

using namespace inchworm;

TEST(Torus, Wrap) {
  EXPECT_CLOSE(wrap(2.0, 1.0), 0.0);
  EXPECT_CLOSE(wrap(1.3, 1.0), 0.3);
  EXPECT_CLOSE(wrap(1.7, 1.0), 0.7);
  EXPECT_CLOSE(wrap(-1.3, 1.0), 0.7);
  EXPECT_CLOSE(wrap(-1.7, 1.0), 0.3);
}

TEST(Torus, CyclicDifference) {
  EXPECT_CLOSE(cyclic_difference(0.3, 0.1, 1.0), 0.2);
  EXPECT_CLOSE(cyclic_difference(0.1, 0.3, 1.0), 0.8);

  EXPECT_CLOSE(cyclic_difference(1.2, 1.1, 1.0), 0.1);
  EXPECT_CLOSE(cyclic_difference(1.1, 1.2, 1.0), 0.9);

  EXPECT_CLOSE(cyclic_difference(2.4, 0.1, 2.0), 0.3);
  EXPECT_CLOSE(cyclic_difference(2.1, 0.4, 2.0), 1.7);

  EXPECT_CLOSE(cyclic_difference(0.5, 2.1, 2.0), 0.4);
  EXPECT_CLOSE(cyclic_difference(0.1, 2.5, 2.0), 1.6);
}

TEST(Torus, CyclicDifferenceConfig) {
  // ------------ tau,   dag, i, bl
  auto d1 = fop_t{0.1, false, 0, 0};
  auto d2 = fop_t{0.8, false, 0, 0};
  auto d3 = fop_t{0.4, false, 0, 1};
  auto d4 = fop_t{0.0, false, 0, 1};

  auto d_dag1 = fop_t{0.3, true, 0, 0};
  auto d_dag2 = fop_t{0.7, true, 0, 0};
  auto d_dag3 = fop_t{0.9, true, 0, 1};
  auto d_dag4 = fop_t{0.5, true, 0, 1};

  config_t config({}, {{"0", 1}, {"1", 1}}, {0.0, 0.5});

  EXPECT_CLOSE(cyclic_difference(d1, config, 1.0), 0.1);
  EXPECT_CLOSE(cyclic_difference(d2, config, 1.0), 0.8);
  EXPECT_CLOSE(cyclic_difference(d3, config, 1.0), 0.9);
  EXPECT_CLOSE(cyclic_difference(d4, config, 1.0), 0.5);

  EXPECT_CLOSE(cyclic_difference(d_dag1, config, 1.0), 0.8);
  EXPECT_CLOSE(cyclic_difference(d_dag2, config, 1.0), 0.2);
  EXPECT_CLOSE(cyclic_difference(d_dag3, config, 1.0), 0.9);
  EXPECT_CLOSE(cyclic_difference(d_dag4, config, 1.0), 0.5);
}

TEST(Torus, CyclicDistance) {
  EXPECT_CLOSE(cyclic_distance(0.1, 0.2, 1.0), 0.1);
  EXPECT_CLOSE(cyclic_distance(0.2, 0.1, 1.0), 0.1);
  EXPECT_CLOSE(cyclic_distance(0.8, 0.4, 1.0), 0.4);
  EXPECT_CLOSE(cyclic_distance(0.4, 0.8, 1.0), 0.4);
  EXPECT_CLOSE(cyclic_distance(0.2, 0.8, 1.0), 0.4);
  EXPECT_CLOSE(cyclic_distance(0.8, 0.2, 1.0), 0.4);
}
