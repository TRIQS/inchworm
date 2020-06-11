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

namespace inchworm::diagram {
  constexpr int verbose = 0;

  // optimization:
  constexpr int smallest_segment   = 4;     // must be 2 or 4, beware.
  constexpr bool remove_xoxo       = true;  // new optimisation 1: IMPORTANT, only works with the option smallest_segment = 4;
  constexpr bool remove_not_finite = true; // new optimisation 2: do not calculate determinant of
                                            // the remainder when one segment is equal to zero (double == 0.0).
                                            // But it might be not as straightforward when doing this comparison for complex values.

  // check if an arch cross a point,
  //
  // example1: cross
  //             p
  //         a______b
  //
  // example2: do not cross
  //     p
  //         a______b
  //
  // example3: do not cross
  //         p
  //         a______b
  //
  // example4: do not cross
  //                p
  //         a______b
  //
  inline bool segment_cross_p(int a, int b, int p) { return (a - p) * (b - p) < 0; }

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
  inline bool arches_cross(int a1, int b1, int a2, int b2) { return (a1 - a2) * (b1 - a2) * (b1 - b2) * (a1 - b2) < 0; }

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
  inline bool arch_crosses_point(int a, int b, int point) { return (a - (point + 0.5)) * (b - (point + 0.5)) < 0.0; }

  // calculate n!
  //
  inline int factorial(int n) {
    if (n > 1)
      return n * factorial(n - 1);
    else
      return 1;
  }

} // namespace inchworm::diagram
