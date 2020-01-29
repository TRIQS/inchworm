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

//#include <inchworm/solver_core.hpp>
#include <inchworm/impurity_product.hpp>
using namespace inchworm;

/*
#include <numeric>
#include <bitset>

#include <triqs/utility/macros.hpp>
#include <triqs/gfs.hpp>
#include <triqs/test_tools/gfs.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/atom_diag/atom_diag.hpp>
#include <inchworm/diagram/diagram.hpp>

using namespace inchworm;
//using propagator_t = block_gf<imtime>;
//using atom_diag = triqs::atom_diag::atom_diag<false>;
using scalar_t = double;
using matrix_t = matrix<scalar_t>;

// Necessary to use atom_diag block structure for the propagator.
triqs::hilbert_space::gf_struct_t find_propagator_struct(triqs::atom_diag::atom_diag<false> const &ad) {
  int n_sub = ad.n_subspaces();
  triqs::hilbert_space::gf_struct_t propagator_struct;

  std::printf("%d: \n", n_sub);
  for (int i = 0; i < n_sub; i++) {
    //int sub_dim = ad.get_subspace_dim(i);
    std::printf("%d ", ad.get_subspace_dim(i));

    std::vector<std::variant<int, std::string>> l(ad.get_subspace_dim(i));
    std::iota(l.begin(), l.end(), 0);
    propagator_struct.push_back(std::make_pair(std::to_string(i), l));
  }
  std::printf("\n\n");

  return propagator_struct;
}


struct propagator_frame {
  std::vector<matrix<dcomplex>> matrices;
  int acc_number;

  propagator_frame(triqs::atom_diag::atom_diag<false> const &ad) : matrices(ad.n_subspaces()), acc_number(0) {
    for (int bl = 0; bl < ad.n_subspaces(); bl++) {
      matrices[bl] = matrix<dcomplex>(ad.get_subspace_dim(bl), ad.get_subspace_dim(bl));
      matrices[bl] = 0;
    }
  }
  propagator_frame &operator+=(propagator_frame U_frame) {
    for (int bl = 0; bl < matrices.size(); bl++) matrices[bl] += U_frame.matrices[bl];
    acc_number++;
    return *this;
  }
  void assign(int bl, matrix<dcomplex> mat) {
    if (acc_number > 1) {
      std::printf("error: cannot assign in an accumalted frame.\n");
      exit(0);
    } else {
      acc_number = 1;
    }
    matrices[bl] += mat;
  }
  void reset() {
    for (int bl = 0; bl < matrices.size(); bl++) matrices[bl] = 0;
  }

  double frobenius_norm() {
    double val = 0;
    for (int bl = 0; bl < matrices.size(); bl++) {
      for (int i = 0; i < first_dim(matrices[bl]); i++) {
        for (int j = 0; j < second_dim(matrices[bl]); j++) {

          double elem = std::abs(matrices[bl](i, j));
          val += elem * elem;
        }
      }
      if (acc_number > 1) {
        std::printf("error: frobenius_norm \n");
        exit(0);
      };
    }
    return std::sqrt(val);
  }

  friend std::ostream &operator<<(std::ostream &out, propagator_frame const &U_frame) {
    out << "propagator_frame (size: " << U_frame.matrices.size() << ")\n";
    for (int bl = 0; bl < U_frame.matrices.size(); bl++) { out << U_frame.matrices[bl] << "\n"; }
    return out;
  }
};

propagator_frame propagator_product(u_tau_t const &U, triqs::atom_diag::atom_diag<false> const &ad, time_diagram_t const &diagram, double tau,
                                    bool use_bare_U) {

  if (not use_bare_U) {
    EXPECTS(U.size() == ad.n_subspaces()); //??? this test does not seems to work????i
  }
  //auto fs = ad.get_fock_states();

  //std::vector<matrix<dcomplex>> U_frame(ad.n_subspaces());
  propagator_frame U_frame(ad);

  for (int initial_bl = 0; initial_bl < ad.n_subspaces(); initial_bl++) {
    int dim = ad.get_subspace_dim(initial_bl);
    //U_frame[initial_bl] = matrix<dcomplex>(dim, dim);
    //U_frame[initial_bl] = 0;
    //std::cout << "U_mat\n" << U_frame[initial_bl] << "\n\n";

    int new_bl = initial_bl;
    for (int i = diagram.size() - 1; i >= 0; i--) {
      //for (auto const op : diagram.op_list) {
      auto op = diagram.op_list[i];
      new_bl  = (op.dag ? ad.cdag_connection(op.linear_index, new_bl) : ad.c_connection(op.linear_index, new_bl));
      //std::printf(" %d  %d   % 4.5f   ", op.linear_index, op.dag, op.tau);
      //std::printf("  old: %d, new: %d \n", initial_bl, new_bl);
      if (new_bl == -1)
        break; /// !!!!!!!!!! important because connection(*, -1) is not correct (should give -1). Need additional optimization. does not take into account consecutive c_i c_i, or cdag_i cdag_i
    }
    //std::printf("bloc: %d, goes to: %d \n", initial_bl, new_bl);

    //for (int i = 0; i < diagram.size(); i++) {
    matrix<dcomplex> new_mat = matrix<dcomplex>(dim, dim);
    new_mat                  = 0;
    if (new_bl != -1) {
      new_bl = initial_bl;

      double dtau = tau - diagram.max_tau();
      //std::cout << "dtau:" << dtau << "\n"
      //          << "tau:" << tau << "\n"
      //          << "diagram.max_tau():" << diagram.max_tau() << "\n";
      //std::cout << "before\n" << U[initial_bl](dtau) << "\n\n";
      //matrix<dcomplex> new_mat = U[initial_bl][0];
      if (use_bare_U) {
        for (int j = 0; j < dim; j++)
          new_mat(j, j) = std::exp(-dtau * ad.get_eigenvalue(initial_bl, j)); // Create time-evolution matrix e^-H(tau-tau_max)
      } else {
        new_mat = U[initial_bl](dtau);
      }
      //std::cout << "after\n" << new_mat << "\n\n";

      for (int i = diagram.size() - 1; i >= 0; i--) {
        // for (auto const op : diagram.op_list) {
        auto op = diagram.op_list[i];
        //std::cout << "just_before\n" << new_mat << "\n\n";
        //std::cout << "new_bl " << new_bl << "\n\n";
        //std::cout << "c? " << ad.c_matrix(op.linear_index, new_bl) << "\n\n";
        //std::cout << "cdag? " << ad.cdag_matrix(op.linear_index, new_bl) << "\n\n";
        //std::cout << "op.dag " << (op.dag ? "true " : "false ") << "\n\n";
        new_mat = (op.dag ? ad.cdag_matrix(op.linear_index, new_bl) * new_mat : ad.c_matrix(op.linear_index, new_bl) * new_mat);
        new_bl  = (op.dag ? ad.cdag_connection(op.linear_index, new_bl) : ad.c_connection(op.linear_index, new_bl));

        if (i < diagram.size() - 1) {
          dtau = op.tau - diagram.op_list[i - 1].tau;
        } else {
          dtau = op.tau;
        }
        //std::cout << "mat:" << matrix<dcomplex>{U[new_bl](dtau)} << "\n";
        //std::cout << "mat:" << new_mat << "\n";
        if (use_bare_U) {
          auto _ = arrays::range();
          for (int j = 0; j < dim; j++) new_mat(_, j) *= std::exp(-dtau * ad.get_eigenvalue(initial_bl, j)); // Time-evolution
        } else {
          new_mat = matrix<dcomplex>{U[new_bl](dtau)} * new_mat;
        }
        std::cout << "mat:" << new_mat << "\n";
      }
      U_frame.assign(new_bl, new_mat);
    }
    std::cout << "mat:" << new_mat << "\n";
    std::printf("\n\n");
  }
  return U_frame;
}

propagator_frame propagator_product(triqs::atom_diag::atom_diag<false> const &ad, time_diagram_t const &diagram, double tau, bool use_bare_U) {
  u_tau_t U;
  return propagator_product(U, ad, diagram, tau, use_bare_U);
}
*/

void print_bin(int v) { std::cout << std::bitset<4>(v); }

void print_ad(triqs::atom_diag::atom_diag<false> const &ad) {
  printf("\n%d \n", ad.get_full_hilbert_space_dim());
  for (auto fock_sp : ad.get_fock_states()) {
    for (auto l : fock_sp) {
      print_bin(l);
      std::printf(" ");
    }
    std::printf("\n");
  }
  std::printf("\n\n");
}

TEST(inchworm, matrix_product) {

  int n_times = 4;
  double beta = 1.0;
  double mu   = 0.0;
  double U    = 8.0;
  double t    = 1.0;
  fundamental_operator_set fops;
  int n_site = 2;

  auto qn_vector = std::vector<triqs::operators::many_body_operator_generic<scalar_t>>();
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
  auto propagator        = u_tau_t{{beta, Fermion, n_times}, propagator_struct};
  //for (int i = 0; i < ad.n_subspaces(); i++) { std::cout << propagator[i] << "\n"; }
  //print_block_gf_first_time(propagator);

  //print_ad(ad);
  //print_block_gf_first_time(propagator);
  std::vector<time_and_indices_t> c    = {{0.00, 0}, {0.05, 0}};
  std::vector<time_and_indices_t> cdag = {{0.01, 0}, {0.08, 0}};
  std::vector<double> split_times      = {0.99};
  time_diagram_t diagram(c, cdag, split_times);

  auto U_frame = propagator_product(ad, diagram, 0.1);
  std::cout << U_frame;
  std::cout << U_frame.frobenius_norm() << "\n";

  assign_identity_to_propagator(propagator,0);
  assign_frame_to_propagator(propagator, U_frame, 1);
  std::cout << U_frame;
  std::cout << U_frame.frobenius_norm() << "\n";

  U_frame = propagator_product(propagator, ad, diagram, 0.2);
  assign_frame_to_propagator(propagator, U_frame, 2);
  std::cout << U_frame;
  std::cout << U_frame.frobenius_norm() << "\n";

  U_frame = propagator_product(propagator, ad, diagram, 0.3);
  assign_frame_to_propagator(propagator, U_frame, 3);
  std::cout << U_frame;
  std::cout << U_frame.frobenius_norm() << "\n";

/* 
  //propagator_product(8.0, ad, diagram, 0.4); //////////ATTENTION CA MARCHE (POURQUOI?)
*/  //print_block_gf_first_time(propagator);
}

MAKE_MAIN
