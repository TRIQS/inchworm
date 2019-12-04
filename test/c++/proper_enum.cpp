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

//#include <triqs/gfs.hpp>
//#include <triqs/h5.hpp>
#include <triqs/test_tools/gfs.hpp>
#include <inchworm/diagram/proper_enum.cpp>

TEST(inchworm, dummy) { 
  std::string diagram ="oxoxxoxxoooxxoxo";
  int split_point = 3;
  int verbose = 0;
  int N_proper = find_proper_diagrams(diagram, split_point, verbose);

  EXPECT_EQ(N_proper, 10832);
}

MAKE_MAIN
