#include "./impurity_product.hpp"

namespace inchworm {

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

  u_frame_t init_propagator_frame(triqs::atom_diag::atom_diag<false> const &ad) {
    u_frame_t u_frame{ad.n_subspaces()};

    for (int bl = 0; bl < ad.n_subspaces(); bl++) {
      u_frame(bl) = matrix<dcomplex>(ad.get_subspace_dim(bl), ad.get_subspace_dim(bl));
      u_frame(bl) = 0;
    }
  }

  // Calculate the Frobenius norm of the u_frame block diagonal matrix:
  double frobenius_norm(u_frame_t u_frame) {
    double val = 0;
    for (int bl = 0; bl < u_frame.size(); bl++) {
      for (int i = 0; i < first_dim(u_frame(bl)); i++) {
        for (int j = 0; j < second_dim(u_frame(bl)); j++) {

          double elem = std::abs(u_frame(bl)(i, j));
          val += elem * elem;
        }
      }
    }
    return std::sqrt(val);
  }


/*
  // Constructor
  // Function to add them, and accumulate.
  // u_frame_t::u_frame_t &operator+=(u_frame_t u_frame) {
  //  for (int bl = 0; bl < matrices.size(); bl++) matrices[bl] += u_frame.matrices[bl];
  //  acc_number++;
  //  return *this;
  //}

  void u_frame_t::assign(int bl, matrix<dcomplex> mat) {
    if (acc_number > 1) {
      std::printf("error: cannot assign in an accumalted frame.\n");
      exit(0);
    } else {
      acc_number = 1;
    }
    matrices[bl] = mat;
  }

  // Set the values to zero
  void u_frame_t::reset() {
    for (int bl = 0; bl < matrices.size(); bl++) matrices[bl] = 0;
    acc_number = 0;
  }

  // Calculate the Frobenius norm of the matrix
  double u_frame_t::frobenius_norm() {
    double val = 0;
    for (int bl = 0; bl < matrices.size(); bl++) {
      for (int i = 0; i < first_dim(matrices[bl]); i++) {
        for (int j = 0; j < second_dim(matrices[bl]); j++) {

          double elem = std::abs(matrices[bl](i, j));
          val += elem * elem;
        }
      }
      if (acc_number > 1) {
        std::printf("error: frobenius_norm...\n");
        exit(0);
      };
    }
    return std::sqrt(val);
  }

  // Printing function.
  std::ostream &operator<<(std::ostream &out, u_frame_t const &u_frame) {
    out << "u_frame_t (size: " << u_frame.matrices.size() << ")\n";
    for (int bl = 0; bl < u_frame.matrices.size(); bl++) { out << u_frame.matrices[bl] << "\n"; }
    return out;
  }
*/

  void init_propagator(u_tau_t &u_tau) { // this assign identity to the first frame (or time) of the propagator.
    for (int bl = 0; bl < u_tau.size(); bl++) u_tau[bl][0] = make_unit_matrix<dcomplex>(first_dim(u_tau[bl][0]));
    return;
  }

  /// Function that calculate the product: u_frame = U(tau_0) op U(tau_1-tau_0) op U(tau_2-tau_1) op U(tau_3-tau_2) ... op U(tau-tau_n)
  /// where op is either c_dag or c operator, depending on the configuration
  /** 
   * @param U Full propagator calculated up until this point.
   * @param ad atom_diag of the system considered here.
   * @param diagram Configuration of the n operators (op) of the present Monte Carlo step.
   * @param tau Time of the u_frame_t calculated here. tau must be greater than any times
   * @param use_bare_U If true, calculate the same product using only the bare propagators. 
   * @return u_frame_t, at time tau, resulting from this product.
   */
  u_frame_t propagator_product(u_tau_t const &U, triqs::atom_diag::atom_diag<false> const &ad, time_diagram_t const &diagram, double tau,
                                      bool use_bare_U) {

    if (not use_bare_U) {
      //EXPECTS(false);
      EXPECTS(U.size() == ad.n_subspaces()); //??? this test does not seems to work????
    }
    //auto fs = ad.get_fock_states();

    u_frame_t u_frame = init_propagator_frame(ad);

    for (int initial_bl = 0; initial_bl < ad.n_subspaces(); initial_bl++) {
      int dim                  = ad.get_subspace_dim(initial_bl);
      int new_bl               = initial_bl;
      matrix<dcomplex> new_mat = matrix<dcomplex>(dim, dim);
      new_mat                  = 0;

      //std::printf("bl %d \n", initial_bl);
      if (diagram.size() == 0) { // if order is zero, use bare propagator and skip the rest of the function
        for (int j = 0; j < dim; j++) {
          new_mat(j, j) = std::exp(-(0-tau) * ad.get_eigenvalue(initial_bl, j));
        }
        u_frame(new_bl) = new_mat;
        continue;
      }

      for (int i = diagram.size() - 1; i >= 0; i--) {
        //std::printf("i %d \n", i);
        auto op = diagram.op_list[i];
        new_bl  = (op.dag ? ad.cdag_connection(op.linear_index, new_bl) : ad.c_connection(op.linear_index, new_bl));
        //std::printf(" %d  %d   % 4.5f   ", op.linear_index, op.dag, op.tau);
        //std::printf("  old: %d, new: %d \n", initial_bl, new_bl);
        if (new_bl == -1)
          break; /// !!!!!!!!!! important because connection(*, -1) is not correct (should give -1). Need additional optimization. does not take into account consecutive c_i c_i, or cdag_i cdag_i
      }
      //std::printf("bloc: %d, goes to: %d \n", initial_bl, new_bl);

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
            auto _ = triqs::arrays::range();
            for (int j = 0; j < dim; j++) new_mat(_, j) *= std::exp(-dtau * ad.get_eigenvalue(initial_bl, j)); // Time-evolution
          } else {
            new_mat = matrix<dcomplex>{U[new_bl](dtau)} * new_mat;
          }
          //std::cout << "mat:" << new_mat << "\n";
        }
        u_frame(new_bl) = new_mat;
      }
      //std::cout << "mat3:" << new_mat << "\n";
      //std::printf("\n\n");
      //std::cout << "u_frame\n" << u_frame << "\n\n";
    }
    //std::cout << "u_frame_end\n" << u_frame << "\n\n";
    return u_frame;
  } // namespace inchworm

  /// If the user do not provide the propagator, use bare propagator instead.
  u_frame_t propagator_product(triqs::atom_diag::atom_diag<false> const &ad, time_diagram_t const &diagram, double tau) {
    u_tau_t U;
    return propagator_product(U, ad, diagram, tau, true);
  }

} // namespace inchworm
