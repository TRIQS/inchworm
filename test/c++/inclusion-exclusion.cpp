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
#include <inchworm/diagram/inclusion-exclusion.cpp>

TEST(inchworm, subdeterminant) { 
  std::vector<time_and_orbital_t> tau1={{1.0,0},{2.0,0},{3.0,0},{4.5,0}};
  std::vector<time_and_orbital_t> tau2={{1.2,0},{1.8,0},{4.6,0},{4.7,0}};
  op_list_t diagram(tau1, tau2);

  hybridization_matrix hyb_mat(diagram);
  
  //printf("det=% 4.8f\n",hyb_mat.extract_det());
  //double value = hyb_mat.extract_det({0,2});
  //printf("det=% 4.8f\n",value);
  EXPECT_CLOSE(1.7968,hyb_mat.extract_det({0,2}));
  EXPECT_CLOSE(-13.43255200,hyb_mat.extract_det({0,2,5,6}));
  EXPECT_CLOSE(0.13697019,hyb_mat.det());
}


TEST(inchworm, segment) { 
  int verbose = 1;
  std::vector<time_and_orbital_t> tau1={{1.0,0},{2.0,0},{3.0,0},{4.5,0}};
  std::vector<time_and_orbital_t> tau2={{1.2,0},{1.8,0},{4.6,0},{4.7,0}};
  op_list_t diagram(tau1, tau2);
  int split_point = 4;
  
  printDiag(split_point, diagram);
  std::vector<k_connected_segment_t> segment_list = determine_segments(split_point, diagram, 1);
  printf("\nsegment_list = %lu\n\n",segment_list.size());
  
  std::vector<combination_of_segments_t> combination_disjoint_list = combine_segments(segment_list, diagram.k_order(), split_point, true);    
  std::vector<combination_of_segments_t> combination_adjacent_list = combine_segments(segment_list, diagram.k_order(), split_point, false);    
  
  if(verbose>0){
    printf("\n\ncombination of disjoint segments:\n\n");
    printDiag(split_point, diagram);  
    for(auto comb : combination_disjoint_list){
      comb.print(segment_list);
      printf("\n");
    }
    printf("\n\ncombination of adjacent segments:\n\n");  
    printDiag(split_point, diagram);
    for(auto comb : combination_adjacent_list){
      comb.print(segment_list);
      printf("\n");
    }
  }
  
  hybridization_matrix hyb_mat(diagram);
  for(int length=2; length<=2*diagram.k_order(); length+=2){
    if(verbose>0) printf("\n############\nsegment length = %d\n",length);
    for(auto seg : segment_list){
      if(seg.size == length){
        calculate_segment(seg.numero, segment_list, combination_disjoint_list, combination_adjacent_list, hyb_mat, diagram.k_order(), 1);
      }
    }
  }
  printf("segment_list.back().value % 4.7f\n",segment_list.back().value);
  printf("segment_list.back().value_without_cuts % 4.7f\n",segment_list.back().value_without_cuts);
}

MAKE_MAIN

