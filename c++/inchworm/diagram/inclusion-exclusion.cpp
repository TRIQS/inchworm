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
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <numeric>

#include <triqs/gfs.hpp>
#include <triqs/utility/macros.hpp>
#include <triqs/utility/itertools.hpp>

#define SMALLEST_SEGMENT 4

struct time_and_orbital_t {
  double tau = 0.;
  int orb    = 0;
};

bool operator<(time_and_orbital_t const &t1, time_and_orbital_t const &t2) { return (t1.tau < t2.tau); }

class diagram_t {

  struct op_t {
    double tau      = 0.;
    bool dag        = false;
    int orb         = 0;
    int order_index = 0;
  };

  public:
  std::vector<op_t> list;
  std::vector<time_and_orbital_t> c_list, cdag_list;
  //std::vector<double> split_points;
  std::vector<double> split_points;
  int k_order() const { return c_list.size(); }

  diagram_t(std::vector<time_and_orbital_t> const &c, std::vector<time_and_orbital_t> const &cdag, std::vector<double> const &split_times)
     : list(2 * c.size()), c_list{c}, cdag_list{cdag} {

    EXPECTS(std::is_sorted(c.begin(), c.end()));
    EXPECTS(std::is_sorted(cdag.begin(), cdag.end()));
    EXPECTS(c.size() == cdag.size());

    split_points.reserve(split_times.size());

    int k_order = c.size();

    for (int i = 0, j = k_order; i < k_order; i++, j++) {
      list[i].tau         = c[i].tau;
      list[i].orb         = c[i].orb;
      list[i].dag         = false;
      list[i].order_index = i;

      list[j].tau         = cdag[i].tau;
      list[j].orb         = cdag[i].orb;
      list[j].dag         = true;
      list[j].order_index = i;
    }

    std::sort(list.begin(), list.end(), [](auto const &x, auto const &y) { return x.tau < y.tau; });

    printf("split points:\n");
    for (auto s_time : split_times) {
      int i;
      for (i = 0; i < list.size(); i++) {
        if (s_time < list[i].tau) break;
      }
      split_points.push_back(i);
      printf("%d  % 4.3f\n", i, s_time);
    }
  }

  std::string string() const {
    std::string diagram_order = "";
    for (auto l : list) {
      if (l.dag)
        diagram_order += "x";
      else
        diagram_order += "o";
    }
    return diagram_order;
  }
};

void print_vector(std::vector<int> const &v) {
  for (auto l : v) { printf("%d ", l); }
  printf("\n");
}

using hybridization_scalar_t = double;
//using hybridization_function_t = triqs::gfs::gf<triqs::gfs::imtime,triqs::gfs::matrix_real_valued>;
using triqs::utility::enumerate;

hybridization_scalar_t hyb_function(hybridization_scalar_t dtau) { // to be changed in the future
  return (2.2 + dtau + 0.7 * dtau * dtau + 0.1 * dtau * dtau * dtau);
}

class hybridization_matrix {

  public:
  using matrix_t = triqs::arrays::matrix<hybridization_scalar_t>;

  matrix_t mat;
  diagram_t diagram;

  //hybridization_matrix(diagram_t const & diagram, hybridization_function_t const & hyb): mat(diagram.k_order(), diagram.k_order()), diagram{diagram}{

  hybridization_matrix(diagram_t const &diagram) : mat(diagram.k_order(), diagram.k_order()), diagram{diagram} {
    //int N = diagram.k_order();

    for (auto [i, c] : enumerate(diagram.c_list))
      for (auto [j, cdag] : enumerate(diagram.cdag_list)) {
        //for(auto const & cdag : diagram.cdag_list){

        double dtau = cdag.tau - c.tau;
        mat(i, j)   = hyb_function(dtau);
        if (SMALLEST_SEGMENT == 4) {
          for (int k = 0; k < diagram.list.size() - 1; k++) {

            if (std::any_of(begin(diagram.split_points), end(diagram.split_points), [k](int l) { return l == k+1; })) continue;

            if (not diagram.list[k].dag and diagram.list[k + 1].dag) { // segment on length 2
              int it          = diagram.list[k + 1].order_index;       /// ERROR: it seems that the definition of c and cdagger are inverted?
              int it_dag      = diagram.list[k].order_index;
              mat(it, it_dag) = 0.;
            } else if (diagram.list[k].dag and not diagram.list[k + 1].dag) { // segment on length 2
              int it_dag      = diagram.list[k + 1].order_index;
              int it          = diagram.list[k].order_index;
              mat(it, it_dag) = 0.;
            }
          }
        }
        //put this below in hyb_funciton at some point:
        //        if(dtau>=0) mat(i,j) = hyb( dtau )( cdag.orb, c.orb );
        //        else mat(i,j) = -hyb( hyb.mesh().domain().beta + dtau )( cdag.orb, c.orb );
      }
  };

  hybridization_scalar_t det() { return determinant(mat); }

  hybridization_scalar_t extract_det(std::vector<int> const &list_of_indices) const {

    int N = list_of_indices.size() / 2;
    EXPECTS(list_of_indices.size() % 2 == 0);
    matrix_t m(N, N);

    std::vector<int> list_of_c(N), list_of_cdag(N); // to refactor

    int i = 0, j = 0;
    for (auto idx : list_of_indices) {
      if (diagram.list[idx].dag)
        list_of_cdag[i++] = diagram.list[idx].order_index;
      else
        list_of_c[j++] = diagram.list[idx].order_index;
    }

    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) { m(i, j) = mat(list_of_cdag[i], list_of_c[j]); }
    }

    //print_vector(list_of_indices);
    //TRIQS_PRINT(m);
    //TRIQS_PRINT(determinant(m));

    return determinant(m);
  }
};

// print diagram and its split point above.
//
void printDiag(diagram_t const &diagram) {
  for (int j = 0; j < diagram.list.size(); j++) {
    for (auto split_point : diagram.split_points) {
      if (split_point - 1 == j)
        printf(" |");
      else
        printf("  ");
    }
  }
  printf("\n");
  for (int j = 0; j < diagram.list.size(); j++) {
    if (diagram.list[j].dag)
      printf("x");
    else
      printf("o");
    if (j < diagram.list.size() - 1) printf("-");
  }
  printf("\n");
}

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
bool segment_cross_p(int a, int b, int p) { return (a - p) * (b - p) < 0.0; }

void printLine(std::vector<int> const &segments_vector) {
  int current_segment = 0;
  std::string string1 = "";
  std::string chars   = "  ";
  for (int i = 0; i < segments_vector.size(); i++) {
    if ((current_segment != segments_vector[i]) or (i == 0)) {
      current_segment = segments_vector[i];
      if (current_segment == -1)
        chars = "  ";
      else if (current_segment == 0)
        chars = "--";
      else
        chars = "==";
      if (string1.size() > 0) {
        string1.pop_back();
        string1 += ' ';
      }
    }
    string1 += chars;
  }
  string1.pop_back();
  printf("%s", string1.c_str());
}

struct k_connected_segment_t {

  int pos1                                  = 0;
  int pos2                                  = 0; // note: by definition here, segment goes from index pos1 to pos2-1
  int size                                  = 0;
  hybridization_scalar_t value              = 0.;
  hybridization_scalar_t value_without_cuts = 0.;

  bool calculated = false;
  int numero      = 0;

  k_connected_segment_t(int p1, int p2, int n) : pos1{p1}, pos2{p2}, numero{n} {
    EXPECTS(pos2 > pos1);
    size = pos2 - pos1;
  }

  void print(int N) {
    std::vector<int> v(N, 0);
    for (int k = pos1; k < pos2; k++) v[k] = 1;
    printLine(v);
  };
};

struct combination_of_segments_t {

  int pos1;
  int pos2; // note: by definition here, the combination of segments goes from index pos1 to pos2-1
  int size;
  bool disjoint = true; // we define disjoint when 2 segments does not touch (by convention, we choose a segment alone to be disjoint too)
  bool adjacent = true; // we define adjacent when all segments touches. If one does not, it is false.
  std::vector<int> list;
  int k_order = 0;
  //int split_point = 0;

  combination_of_segments_t(k_connected_segment_t const &seg0, diagram_t const &diagram) : pos1{seg0.pos1}, pos2{seg0.pos2}, size{seg0.size} {
    list.reserve(
       diagram
          .k_order()); //If smallest segment is length 2, we know that this is the maximum number of segments in a combination. If the smallest is 4, then it becomes k_order/2.
    list.push_back(seg0.numero);
  };

  void append(k_connected_segment_t const &seg1, diagram_t const &diagram) {
    if (seg1.pos1 != pos2)
      adjacent = false;
    else if (std::any_of(begin(diagram.split_points), end(diagram.split_points), [j = pos2](int i) { return i == j; }))
      adjacent =
         false; // if the previous combination of segments already ends at a split point, adding another one will make this combination not adjacent anymore.
    else
      disjoint = false; // note, we consider that even if two segments touch at the split point, the combination is still disjoint.

    size = seg1.pos2 - pos1;
    pos2 = seg1.pos2;
    list.push_back(seg1.numero);
  };

  void print(std::vector<k_connected_segment_t> const &segment_list, diagram_t const &diagram) {
    std::vector<int> num_vector(diagram.list.size(), 0);

    for (int j = 0; j < list.size(); j++) {
      //printf("list[j] = %d, size = %d\n", list[j], list.size());
      //printf("list[j].pos1 = %d,  list[j].pos2 = %d\n", segment_list[list[j]].pos1, segment_list[list[j]].pos2);
      for (int k = segment_list[list[j]].pos1; k < segment_list[list[j]].pos2; k++) num_vector[k] = j + 1;
    }
    printLine(num_vector);
    //std::cout << num_vector << "\n";
    //for(auto v:num_vector) printf("v = %d\n", v);
  };
};

std::vector<k_connected_segment_t> determine_segments(diagram_t const &diagram, int verbose = 0) {

  std::vector<k_connected_segment_t> list;
  int N = diagram.list.size();

  //EXPECTS(split_point >= 0);
  //EXPECTS(split_point < N);

  int N_segment = 0;
  for (int i = 0; i < N - 1; i++)                           //starting position of segment
    for (int a = SMALLEST_SEGMENT; a < N - i + 1; a += 2) { //length of segment

      if (std::any_of(diagram.split_points.begin(), diagram.split_points.end(),
                      [i, a](auto &split_point) { return segment_cross_p(i, i + a, split_point); }))
        continue;
      //if(segment_cross_p(i, i+a, split_point)) {//continue;
      //printf("#### %d %d %d ####\n", i,i+a);
      int Ndag = 0;
      for (int j = i; j < i + a; j++)
        if (diagram.list[j].dag) Ndag++;
      if (2 * Ndag == a) // check if same number of cdag an c in the segment starting at i and ending before i+a
      {
        //printf("numero = %d \n", N_segment);
        k_connected_segment_t seg(i, i + a, N_segment++);
        list.push_back(seg);
        if (verbose > 0) {
          seg.print(N);
          //if((a==4) and (diagram.list[i].dag == diagram.list[i+2].dag ) ) printf("  will always be zero, to be ignore in future");
          //printf("\nnumero=%d\n\n",seg.numero);
          printf("\n");
        }
      }
    }

  //lastly, put the last segment (this one is k-connected and not fully connected. So we bypass the condition that it should not cross the split point:
  k_connected_segment_t seg(0, N, N_segment++);
  list.push_back(seg);
  seg.print(N);
  return list;
}

//std::pair<std::vector<combination_of_segments_t>,std::vector<combination_of_segments_t>>
std::vector<combination_of_segments_t> combine_segments(std::vector<k_connected_segment_t> const &segment_list, diagram_t const &diagram,
                                                        bool search_disjoint) {

  int n_seg = segment_list.size();

  std::vector<int> start_index_list = {0, 0};
  std::vector<combination_of_segments_t> combination_list; // return value (pair[1])

  for (int j = 0; j < n_seg; j++) { combination_list.push_back(combination_of_segments_t(segment_list[j], diagram)); }
  for (int number_of_segment = 2; number_of_segment <= diagram.k_order(); number_of_segment++) {
    start_index_list.push_back(combination_list.size());

    int i1 = start_index_list.size() - 2;
    int i2 = start_index_list.size() - 1;
    for (int prev_index = start_index_list[i1]; prev_index < start_index_list[i2]; prev_index++) {

      combination_of_segments_t previous_combination = combination_list[prev_index];
      //if we do not search for disjoint combination, we search for adjacent combination, only. We do not need the ones that are neither.

      for (auto additional_segment : segment_list) {
        if (additional_segment.pos1 >= previous_combination.pos2) {
          combination_of_segments_t new_combination = previous_combination;
          new_combination.append(additional_segment, diagram);

          if (search_disjoint) {
            if (not new_combination.disjoint) continue;
          } else {
            if (not new_combination.adjacent) continue;
          } //important brackets

          combination_list.push_back(new_combination);
        }
      }
    }
  }
  return combination_list;
}

std::vector<int> remove_segment_from_list(std::vector<int> const &list, int segment_min, int segment_size) {
  std::vector<int> output;
  output.reserve(list.size() - segment_size);

  int segment_max = segment_min + segment_size;
  for (auto l : list)
    if ((l < segment_min) or (l >= segment_max)) { output.push_back(l); }

  EXPECTS(output.size() == list.size() - segment_size)
  return output;
}

void calculate_segment(int segment_numero,
                       std::vector<k_connected_segment_t> &segments_list, // not const: modified
                       std::vector<combination_of_segments_t> const &combination_disjoint_list,
                       std::vector<combination_of_segments_t> const &combination_adjacent_list, hybridization_matrix const &hyb_mat, int k_order,
                       int verbose = 0, bool special = false) {

  segments_list[segment_numero].calculated = true;
  if (verbose > 1) { segments_list[segment_numero].print(2 * k_order); }

  std::vector<int> range_of_vertex(segments_list[segment_numero].size);
  std::iota(range_of_vertex.begin(), range_of_vertex.end(), segments_list[segment_numero].pos1);

  segments_list[segment_numero].value += hyb_mat.extract_det(range_of_vertex);

  for (auto subs : combination_disjoint_list) {
    if ((not special) and not((segments_list[segment_numero].pos1 <= subs.pos1) and (segments_list[segment_numero].pos2 > subs.pos2))) continue;
    if (special and (subs.list.size() == 1)
        and ((segments_list[segment_numero].pos1 == subs.pos1) and (segments_list[segment_numero].pos2 == subs.pos2)))
      continue; // this is tricky, might have to change this at some point

    hybridization_scalar_t value = 1.0;
    std::vector<int> range_of_subvertex(range_of_vertex);
    int signe_of_parcollet_charlebois = 1;

    for (auto sub_segment_numero : subs.list) {
      k_connected_segment_t seg = segments_list[sub_segment_numero];
      EXPECTS(seg.calculated);

      value *= -seg.value;

      range_of_subvertex = remove_segment_from_list(range_of_subvertex, seg.pos1, seg.size);

      if (seg.size % 4 != 0)
        if ((seg.pos2 - segments_list[segment_numero].pos1) % 2 == 1) signe_of_parcollet_charlebois *= -1;
    }

    if (range_of_subvertex.size() > 0) {
      hybridization_scalar_t det1 = hyb_mat.extract_det(range_of_subvertex);
      value *= signe_of_parcollet_charlebois * det1;
    }
    segments_list[segment_numero].value += value;
  }

  segments_list[segment_numero].value_without_cuts = segments_list[segment_numero].value;

  for (auto cuts : combination_adjacent_list) {
    if (segments_list[segment_numero].pos1 == cuts.pos1)
      if (segments_list[segment_numero].pos2 == cuts.pos2)
        if (cuts.list.size() > 1) {

          hybridization_scalar_t value = 1.0;

          for (auto sub_segment_numero : cuts.list) {
            k_connected_segment_t seg = segments_list[sub_segment_numero];
            EXPECTS(seg.calculated);

            value *= -seg.value_without_cuts;
          }
          segments_list[segment_numero].value -= value;
        }
  }
  if (verbose > 1) printf("   % 15.8f      % 15.8f\n", segments_list[segment_numero].value, segments_list[segment_numero].value_without_cuts);
}

hybridization_scalar_t inclusion_exclusion(diagram_t const &diagram, int verbose = 0) {

  if (verbose > 1) printDiag(diagram);
  std::vector<k_connected_segment_t> segment_list = determine_segments(diagram, verbose);
  if (verbose) printf("\nsegment number = %lu\n\n", segment_list.size());

  std::vector<combination_of_segments_t> combination_disjoint_list = combine_segments(segment_list, diagram, true);
  std::vector<combination_of_segments_t> combination_adjacent_list = combine_segments(segment_list, diagram, false);

  if (verbose > 1) {
    printf("\n\ncombination of disjoint segments:\n\n");
    printDiag(diagram);
    for (auto comb : combination_disjoint_list) {
      comb.print(segment_list, diagram);
      printf("\n");
    }
    printf("\n\ncombination of adjacent segments:\n\n");
    printDiag(diagram);
    for (auto comb : combination_adjacent_list) {
      comb.print(segment_list, diagram);
      printf("\n");
    }
  }

  hybridization_matrix hyb_mat(diagram);

  for (int length = SMALLEST_SEGMENT; length <= 2 * diagram.k_order(); length += 2) {
    if (verbose > 1) printf("\n############\nsegment length = %d\n", length);
    for (auto seg : segment_list) {
      if (seg.size == length) {
        bool special = false;
        if (length == 2 * diagram.k_order()) special = true;
        calculate_segment(seg.numero, segment_list, combination_disjoint_list, combination_adjacent_list, hyb_mat, diagram.k_order(), verbose,
                          special);
      }
    }
  }

  if (verbose > 0) {
    printf("\n## diagram = '%s'\n", diagram.string().c_str());
    printf("kOrder = %d\n", diagram.k_order());
    printf("number of segments = %lu\n", segment_list.size());
    printf("number of adjacent combinations = %lu\n", combination_adjacent_list.size());
    printf("number of disjoint combinations = %lu\n", combination_disjoint_list.size());
  }

  return segment_list.back().value;
}
