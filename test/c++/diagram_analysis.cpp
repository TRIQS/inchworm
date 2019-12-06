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

/*
TEST(inchworm, dummy) { 
  std::string diagram ="oxoxxoxxoooxxoxo";
  int split_point = 3;
  int verbose = 0;
  int N_proper = find_proper_diagrams(diagram, split_point);

  EXPECT_EQ(N_proper, 10832);
}*/

TEST(inchworm, proper_enum1) { 
  //int verbose = 2;

  //std::vector<double> split_times = {1.9};
  //std::vector<double> tau1={1.0,2.0,3.0,4.5};
  //std::vector<double> tau2={1.2,1.8,4.6,4.7};
  std::vector<double> tau1={0.4344,0.1,0.3,0.5,0.75,0.9};
  std::vector<double> tau2={0,0.3452,0.21,0.45,0.69,0.81};
  std::vector<double> split_times = {0.8};
  
  std::sort (tau1.begin(), tau1.end());
  std::sort (tau2.begin(), tau2.end());
  
  std::vector<time_and_orbital_t> c, cdag;
  for(auto t : tau1) c.push_back({t,0});
  for(auto t : tau2) cdag.push_back({t,0});
  time_diagram_t diagram(c,cdag,split_times);
  
  //int N_proper = find_proper_diagrams(diagram);
  hybridization_scalar_t  value_proper = proper_enum(diagram);
  std::printf("c_k = % 4.6f\n", value_proper);
  
  hybridization_scalar_t  value_inclus = inclusion_exclusion(diagram);
  std::printf("c_k = % 4.6f\n\n\n", value_inclus);
  
  std::printf("proper-enum         c_k = % 4.6f\n",   value_proper);
  std::printf("inclusion-exclusion c_k = % 4.6f\n\n", value_inclus);
    
  //EXPECT_EQ(N_proper, 13);
}








/*
TEST(inchworm, subdeterminant) { 
  std::vector<time_and_orbital_t> tau1={{1.0,0},{2.0,0},{3.0,0},{4.5,0}};
  std::vector<time_and_orbital_t> tau2={{1.2,0},{1.8,0},{4.6,0},{4.7,0}};
  time_diagram_t diagram(tau1, tau2);

  hybridization_matrix hyb_mat(diagram,0);
  
  //std::printf("det=% 4.8f\n",hyb_mat.extract_det());
  //double value = hyb_mat.extract_det({0,2});
  //std::printf("det=% 4.8f\n",value);
  EXPECT_CLOSE(1.7968,hyb_mat.extract_det({0,2}));
  EXPECT_CLOSE(-13.43255200,hyb_mat.extract_det({0,2,5,6}));
  EXPECT_CLOSE(0.13697019,hyb_mat.det());
  EXPECT_CLOSE(hyb_mat.det(),hyb_mat.extract_det({0,1,2,3,4,5,6,7}));
}


TEST(inchworm, segment) { 
  int verbose = 2;
  std::vector<time_and_orbital_t> tau1={{1.0,0},{2.0,0},{4.0,0},{4.5,0}};
  std::vector<time_and_orbital_t> tau2={{0.8,0},{1.8,0},{3.6,0},{4.7,0}};
  std::vector<double> split_times = {2.2};
  time_diagram_t diagram(tau1, tau2, split_times);
  //int split_point = 4;
  
  hybridization_scalar_t  value = inclusion_exclusion(diagram);
  
  std::printf("c_k = % 4.7f\n\n", value);
}

TEST(inchworm, inclusion_exclusion1) { 
  //int verbose = 1;
  //std::vector<double> tau1={0.4344,0.1,0.3,0.5,0.75,0.9,0.55,0.566,0.33,.4959594,.4494929,.12349512,.62343,0.123412444,0.2134444,.99949941,1.1,1.23,1.45,2.3,1.6,4.3,1.222,2.98,3.1244};
  //std::vector<double> tau2={0,0.3452,0.21,0.45,0.69,0.81,0.998,0.122,0.833,0.4934,.210342134,.210343,.02134,.0030404,.02142430,0.1111,1.2,1.3,1.4,3.4,2.34,3.11,2.9,1.99,3.098};

  std::vector<double> tau1={0.4344,0.1,0.3,0.5,0.75,0.9,0.55,0.566,0.33,.4959594,.4494929,.12349512,.62343,0.123412444,0.2134444,.99949941};
  std::vector<double> tau2={0,0.3452,0.21,0.45,0.69,0.81,0.998,0.122,0.833,0.4934,.210342134,.210343,.02134,.0030404,.02142430,0.1111};
  //std::vector<double> split_times = {0.82,0.22,0.32,0.52};
  std::vector<double> split_times = {0.82};
  
  //std::vector<double> tau1={0.4344,0.1,0.3, 0.9};
  //std::vector<double> tau2={0.0,0.3452,0.21,0.7};


  //std::vector<double> tau1={1.0,2.0,3.0,4.5};
  //std::vector<double> tau2={1.23421,1.8,4.6,4.7};
  
  std::sort (tau1.begin(), tau1.end());
  std::sort (tau2.begin(), tau2.end());
  
  std::vector<time_and_orbital_t> c, cdag;
  for(auto t : tau1) c.push_back({t,0});
  for(auto t : tau2) cdag.push_back({t,0});
  time_diagram_t diagram(c,cdag,split_times);
  
  hybridization_scalar_t  value = inclusion_exclusion(diagram);
  std::printf("c_k = % 4.6f\n\n", value);
    
  EXPECT_NEAR(-8458508.82790539,value, 1e-7);
}
*/

MAKE_MAIN

