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
#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <numeric>

#include <triqs/gfs.hpp>
#include <triqs/utility/macros.hpp>
#include <triqs/utility/itertools.hpp>

constexpr int verbose = 3;

// optimization:
constexpr int smallest_segment   = 2;     // must be 2 or 4, beware.
constexpr bool remove_xoxo       = false; // new optimisation 1: IMPORTANT, only works with the option smallest_segment = 4;
constexpr bool remove_not_finite = false; // new optimisation 2: do not calculate determinant of
                                          // the remainder when one segment is equal to zero (double == 0.0).
                                          // But it might be not as straightforward when doing this comparison for complex values.

// check if an arch cross a point,
//
// example1: cross
//             p.
//         a______b
//
// example2: do not cross
//     p.
//         a______b
//
bool segment_cross_p(int a, int b, int p) { return (a - p) * (b - p) < 0; }

// check if two arches cross,
//  i.e. if one end of one arch arrive in the middle of the other arch.
//
// example1: cross
// ____________
//         ________
//
// example2: do not cross
// ___
//         ________
//
// example3: do not cross
//           ___
//         ________
//
bool segment_cross(int a1, int b1, int a2, int b2) { return (a1 - a2) * (b1 - a2) * (b1 - b2) * (a1 - b2) < 0; }

// check if an arch cross a point,
//
// example1: cross
//             |
//         ________
//
// example2: do not cross
//     |
//         ________
//
bool segment_cross_point(int a, int b, int point) { return (a - (point + 0.5)) * (b - (point + 0.5)) < 0.0; }

// calculate n!
//
int factorial(int n) {
  if (n > 1)
    return n * factorial(n - 1);
  else
    return 1;
}
