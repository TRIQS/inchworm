/*******************************************************************************
 *
 * inchworm: A TRIQS based impurity solver
 *
 * Copyright (c) 2019 The Simons foundation
 *   authors: Nils Wentzell
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

#include <numeric>
#include <inchworm/solver_core.hpp>

#include <triqs/gfs.hpp>
//#include <triqs/h5.hpp>
#include <triqs/test_tools/gfs.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/atom_diag/atom_diag.hpp>
//#include <triqs/arrays/blas_lapack/dot.hpp>

#include <inchworm/diagram/diagram.hpp>

using namespace inchworm;
using propagator_t = block_gf<imtime>;
//using atom_diag = triqs::atom_diag::atom_diag<false>;

//block_gf =
triqs::hilbert_space::gf_struct_t find_propagator_struct(triqs::atom_diag::atom_diag<false> const &ad) {
  int n_sub = ad.n_subspaces();
  triqs::hilbert_space::gf_struct_t propagator_struct;

  std::printf("%d: \n", n_sub);
  for (int i = 0; i < n_sub; i++) {
    //int sub_dim = ad.get_subspace_dim(i);
    std::printf("%d ", ad.get_subspace_dim(i));

    std::vector<std::variant<int, std::string>> l(ad.get_subspace_dim(i));
    std::iota(l.begin(), l.end(), 0);
    //for(auto t : l){
    //  std::cout << t << " ";
    //}
    //std::printf("\n\n");
    //auto test = std::make_pair( std::to_string(i), l);
    //std::cout << test.second << " " << test.first << "\n";
    propagator_struct.push_back(std::make_pair(std::to_string(i), l));
  }
  std::printf("\n\n");

  return propagator_struct;
}

void print_block_gf_first_time(block_gf<imtime> const &x) {
  for (int i = 0; i < x.size(); i++) std::cout << x[i][0] << "\n";
  std::printf("\n\n");
}

// todo:
// create an intermediary object for a propagator U for calculation and result (tomorrow).
void propagator_product(propagator_t const &U, triqs::atom_diag::atom_diag<false> const &ad, time_diagram_t const &diagram)
{ 
  for(auto const op: diagram.op_list){
    std::cout << op.dag << "\n";
  }
  std::printf("\n\n");
}

TEST(inchworm, matrix_product) {

  double beta = 10.0;
  double mu   = 0.0;
  double U    = 8.0;
  double t    = 1.0;
  fundamental_operator_set fops;
  int n_site = 2;

  auto qn_vector = std::vector<triqs::operators::many_body_operator_generic<double>>();
  auto h         = 0 * (n("up", 0));
  auto n_tot     = 0 * (n("up", 0));

  for (int i = 0; i < n_site; i++) {
    fops.insert("up", i);
    fops.insert("dn", i);

    h += U * (n("up", i) * (n("dn", i)));
    h -= mu * (n("up", i) + n("dn", i));
    for (int j = 0; j < n_site; j++) {
      if (i != j) h -= t * (c_dag("up", i) * (c("up", j)) + c_dag("dn", i) * (c("dn", j)));
    }
    n_tot += n("up", i) + n("dn", i);
  }
  qn_vector.push_back(n_tot);

  auto ad                = triqs::atom_diag::atom_diag<false>(h, fops, qn_vector);
  auto propagator_struct = find_propagator_struct(ad);

  auto propagator = block_gf<imtime>{{beta, Fermion, 2}, propagator_struct};
  //for (int i = 0; i < ad.n_subspaces(); i++) { std::cout << propagator[i] << "\n"; }
  print_block_gf_first_time(propagator);

  std::vector<time_and_orbital_t> c        = {{0.1,0}, {0.5,1},  {0.6,0}};
  std::vector<time_and_orbital_t> cdag        = {{0.0,0}, {0.11,1}, {0.51,1}};
  std::vector<double> split_times = {0.99};
  time_diagram_t diagram(c, cdag, split_times);
  
  propagator_product(U,ad,diagram);
}

MAKE_MAIN;
