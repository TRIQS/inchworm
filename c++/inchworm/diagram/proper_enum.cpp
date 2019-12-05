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

#include <iostream>
#include <string>
#include <utility>

#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <assert.h> // will be removed later
#include <vector>

// find parity of a permutation by evaluating
// by evaluating the parity of all cycle (or orbits)
//
int find_parity(std::vector<int> const &permutation, int n) {
  int parity = 1;
  //bool visited[n];  //vectorize
  std::vector<bool> visited(n, false);

  for (int ii = 0; ii < n; ii++) {
    if (not visited[ii]) continue;
    int jj = permutation[ii];
    while (jj != ii) {
      parity *= -1; //change for addiiton?
      visited[jj] = true;
      jj          = permutation[jj];
    }
  }
  return parity;
}

// calculate n!
//
int factorial(int n) {
  if (n > 1)
    return n * factorial(n - 1);
  else
    return 1;
}

// return a vector of position where the char "c" appear in the string "str"
//
std::vector<int> find_positions(std::string const &str, char c) {
  std::vector<int> v;
  for (int i = 0; i < str.size(); i++)
    if (str[i] == c) v.push_back(i);
  return v;
}

// print diagram and its split point above.
//
void printDiag(int split_point, std::string const &diagram) {
  for (int j = 0; j < diagram.size(); j++) {
    if (diagram.size() - split_point - 1 == j) {
      printf(" |");
    } else {
      printf("  ");
    }
  }
  printf("\n");
  for (int j = 0; j < diagram.size(); j++) {
    printf("%c", diagram[j]);
    if (j < diagram.size() - 1) printf("-");
  }
  printf("\n");
}

// print the arch from a to b with different character depending on the 'status' of the arch.
//
void printArch(int a, int b, int k_order, char char1 = '.') {
  int ii;
  if (a > b) {
    int tmp = b;
    b       = a;
    a       = tmp;
  } else if (a == b) {
    printf("error a==b");
    exit(1);
  }

  for (ii = 0; ii < 2 * k_order; ii++) {
    if ((ii >= a) and (ii < b))
      printf("%c%c", char1, char1);
    else if ((ii == b))
      printf("%c ", char1);
    else
      printf("  ");
  }
  printf("\n");
}

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

// print graph including diagram and all the arches for
// a specific permutation. Different char are used to represent different status of the arch.
//
void print_graph(const std::vector<int> &permutation, const std::vector<bool> &cross_split_point, std::vector<bool> &visited,
                 const std::vector<int> &pos_o, const std::vector<int> &pos_x, const int k_order, const int split_point, std::string diagram) {

  printDiag(split_point, diagram);

  int arch;
  for (arch = 0; arch < k_order; arch++) { //as much arches than pair of vertices
    char char1 = '.';
    if (cross_split_point[arch]) char1 = '_';
    if (visited[arch]) char1 = '=';
    printArch(pos_o[arch], pos_x[permutation[arch]], k_order, char1);
  }
  printf("\n");
}

// Deep First Search (DFS) algorithm to search for every connected arch.
// This recursive function will call itself until there is no more
// free arch to visit (stored in variable visited).
//
void grow_pile(int arch, std::vector<int> const &permutation, std::vector<bool> &visited, std::vector<int> const &pos_o,
               std::vector<int> const &pos_x, int k_order) {

  //connexion_pile.push_back(arch);
  int o = pos_o[arch];
  int x = pos_x[permutation[arch]];

  for (int new_arch = 0; new_arch < k_order; new_arch++) {
    if (visited[new_arch]) continue;

    if (segment_cross(o, x, pos_o[new_arch], pos_x[permutation[new_arch]])) {
      visited[new_arch] = true;
      grow_pile(new_arch, permutation, visited, pos_o, pos_x, k_order);
    }
  }
  //connexion_pile.pop_back();
  visited[arch] = true;
}

// DFS algorithm to search for every connected arch.
// Everything is set up such that we can call the recursive DFS algorithm
// "grow_pile" function. We need to first make a variable that cross
// that keep track of the arches that cross the split_point, we name
// this vector cross_split_point_pile.
//
bool test_diagram_connection(int split_point, const std::vector<int> &permutation, const std::vector<int> &pos_o, const std::vector<int> &pos_x,
                             int k_order, std::string diagram, int verbose = 0) {

  std::vector<bool> cross_split_point(k_order, false); // for graphic purpose only (a changer)
  std::vector<bool> visited(k_order, false);

  std::vector<int> cross_split_point_pile;
  cross_split_point_pile.reserve(k_order); // we know that the pile will not grow bigger than the number of arch = k_order

  for (int arch = 0; arch < k_order; arch++) {
    int o = pos_o[arch];
    int x = pos_x[permutation[arch]];
    if (segment_cross_point(o, x, 2 * k_order - split_point - 1)) {
      cross_split_point[arch] = true;
      cross_split_point_pile.push_back(arch);
    }
  }

  if (verbose > 1) print_graph(permutation, cross_split_point, visited, pos_o, pos_x, k_order, split_point, diagram);

  for (auto arch : cross_split_point_pile) {
    if (visited[arch]) continue;
    visited[arch] = true;
    grow_pile(arch, permutation, visited, pos_o, pos_x, k_order);
  }

  if (verbose > 1) print_graph(permutation, cross_split_point, visited, pos_o, pos_x, k_order, split_point, diagram);

  return std::all_of(visited.begin(), visited.end(), [](bool v) { return v; });
}

// Look at all the arches combinaitions (represented by a permutation vector)
// for a given diagram definition. Along with the split_point, we can
// use the function "test_diagram_connection" to define if a diagram is
// proper or not.
//
int find_proper_diagrams(std::string const &diagram, int split_point, int verbose = 0) {
  int k_order = diagram.size() / 2;
  std::vector<int> permutation(k_order);
  for (int ii = 0; ii < k_order; ii++) { permutation[ii] = ii; }

  std::vector<int> pos_o = find_positions(diagram, 'o');
  std::vector<int> pos_x = find_positions(diagram, 'x');

  if (pos_o.size() != pos_x.size()) { // ENSURE
    printf("error, must have same o and x in diagram.\n");
    exit(1);
  }

  int NN = 0, N_proper = 0;
  do {
    NN += 1;
    //int parity = find_parity(permutation,k_order);
    if (test_diagram_connection(split_point, permutation, pos_o, pos_x, k_order, diagram, verbose)) { N_proper += 1; }
  } while (std::next_permutation(permutation.begin(), permutation.end()));

  if (verbose) {
    printf("## diagram = '%s'\n", diagram.c_str());
    printf("k_order = %d\n", k_order);
    printf("number of diagram = %d\n", NN);
    printf("number of proper diagram = %d\n", N_proper);
  }
  return N_proper;
}

/*
int main()
{
  //test1();
  //std::string diagram ="oxoxxxooxo";
  //std::string diagram ="oxoxxoxxoooxxoxo";
  std::string diagram = "oxoxoooxooxxxxxoxoxo";
  int split_point = 3;
  int verbose = 0;

  printf("%d! =% d permutations\n",(int) diagram.size()/2,factorial((int) diagram.size()/2)) ;
  find_proper_diagrams(diagram, split_point, verbose);


  return 0;
}
*/
