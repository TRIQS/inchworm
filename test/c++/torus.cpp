/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2020 Simons Foundation
 *   author: N. Wentzell
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include <nda/gtest_tools.hpp>

#include <inchworm/torus.hpp>
#include <inchworm/config.hpp>

using namespace inchworm;

TEST(Torus, Wrap) {
  EXPECT_CLOSE(wrap( 2.0, 1.0), 0.0);
  EXPECT_CLOSE(wrap( 1.3, 1.0), 0.3);
  EXPECT_CLOSE(wrap( 1.7, 1.0), 0.7);
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

  auto d_dag1 = fop_t{0.4, true, 0, 0};
  auto d_dag2 = fop_t{0.7, true, 0, 0};
  auto d_dag3 = fop_t{0.8, true, 0, 1};

  config_t config({}, {{"0", 1}, {"1", 1}}, {0.0});
  config.d_bl_list = {{d1, d2},{d3}};
  config.d_dag_bl_list = {{d_dag1, d_dag2},{d_dag3}};

  EXPECT_CLOSE(cyclic_difference(d1, config, 1.0), 0.7);
  EXPECT_CLOSE(cyclic_difference(d2, config, 1.0), 0.1);
  EXPECT_CLOSE(cyclic_difference(d3, config, 1.0), 0.6);

  EXPECT_CLOSE(cyclic_difference(d_dag1, config, 1.0), 0.3);
  EXPECT_CLOSE(cyclic_difference(d_dag2, config, 1.0), 0.9);
  EXPECT_CLOSE(cyclic_difference(d_dag3, config, 1.0), 0.4);
}

TEST(Torus, CyclicDistance) {
  EXPECT_CLOSE(cyclic_distance(0.1, 0.2, 1.0), 0.1);
  EXPECT_CLOSE(cyclic_distance(0.2, 0.1, 1.0), 0.1);
}
