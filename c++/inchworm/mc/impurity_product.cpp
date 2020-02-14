#include "./impurity_product.hpp"

namespace inchworm {
  //
  triqs::hilbert_space::gf_struct_t find_propagator_struct(atom_diag const &ad) {
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

  //
  u_frame_t make_zero_propagator_frame(atom_diag const &ad) {
    u_frame_t u_frame(ad.n_subspaces());

    for (int bl = 0; bl < ad.n_subspaces(); bl++) {
      u_frame[bl] = matrix_t(ad.get_subspace_dim(bl), ad.get_subspace_dim(bl)); //  use zeros<> ?? check
      u_frame[bl] = 0;
    }
    return u_frame;
  }

  //
  u_frame_t make_bare_propagator_frame(atom_diag const &ad, double tau) {
    u_frame_t u_frame(ad.n_subspaces());

    for (int bl = 0; bl < ad.n_subspaces(); bl++) {
      int dim     = ad.get_subspace_dim(bl);
      u_frame[bl] = matrix_t(dim, dim); //  use zeros<> ?? check
      u_frame[bl] = 0;
      for (int j = 0; j < dim; j++) u_frame[bl](j, j) = std::exp(-tau * ad.get_eigenvalue(bl, j));
    }
    return u_frame;
  }

  // Calculate the Frobenius norm of the u_frame block diagonal matrix:
  double frobenius_norm(u_frame_t const &u_frame) {
    double val = 0;
    for (auto const &B : u_frame) {
      double norm = frobenius_norm(B);
      val += norm * norm;
      //double norm = trace(B);
      //val += norm ;
    }
    return std::sqrt(val);
    //return val;
  }

  // mettre dans une classe
  u_tau_t make_propagator(atom_diag const &h_diag, int n_tau) {
    // this assign identity to the first frame (or time) of the propagator.
    // Build the propagator
    auto propagator_struct = find_propagator_struct(h_diag);
    auto u_tau             = u_tau_t{{1, Fermion, n_tau}, propagator_struct};
    //auto u_frame = u_frame_t{h_diag};

    for (auto &block : u_tau) block[0] = 1; //make_unit_matrix<scalar_t>(first_dim(u_tau[bl][0])); // for auto
    return u_tau;
  }

  //
  u_frame_t propagator_product(atom_diag const &ad, time_diagram_t const &diagram, double tau, u_tau_t const *const u_tau_p) {

    if (diagram.size() == 0) return make_bare_propagator_frame(ad, tau);

    u_frame_t u_frame = make_zero_propagator_frame(ad);
    for (int initial_bl = 0; initial_bl < ad.n_subspaces(); initial_bl++) {
      int dim      = ad.get_subspace_dim(initial_bl);
      int new_bl   = initial_bl;
      auto new_mat = matrix_t{};

      for (int i = 0; (i < diagram.size()) and (new_bl != -1); i++) {
        auto const &op = diagram.op_list[i];
        new_bl         = (op.dag ? ad.cdag_connection(op.linear_index, new_bl) : ad.c_connection(op.linear_index, new_bl));
      }
      //std::printf("bloc: %d, goes to: %d \n", initial_bl, new_bl);

      if (new_bl == -1) continue;
      new_bl = initial_bl;

      double dtau = diagram.min_tau();
      if (u_tau_p)
        new_mat = (*u_tau_p)[initial_bl](dtau);
      else {
        new_mat = matrix_t(dim, dim); //zeros?
        new_mat = 0;
        for (int j = 0; j < dim; j++)
          new_mat(j, j) = std::exp(-dtau * ad.get_eigenvalue(initial_bl, j)); // Create time-evolution matrix e^-H(tau-tau_max)
      }
      //std::cout << "after\n" << new_mat << "\n\n";

      for (int i = 0; i < diagram.size(); i++) {
        auto const &op = diagram.op_list[i];
        //std::cout << "new_mat 1: " << new_bl << " \n" << new_mat << "\n\n";
        //std::cout << (op.dag ? ad.cdag_matrix(op.linear_index, new_bl) : ad.c_matrix(op.linear_index, new_bl)) << "\n\n";
        new_mat = (op.dag ? ad.cdag_matrix(op.linear_index, new_bl) : ad.c_matrix(op.linear_index, new_bl)) * new_mat;
        new_bl  = (op.dag ? ad.cdag_connection(op.linear_index, new_bl) : ad.c_connection(op.linear_index, new_bl));

        //std::cout << "new_mat 2: " << new_bl << " \n" << new_mat << "\n\n";
        dtau = (i == (diagram.size() - 1) ? tau : diagram.op_list[i + 1].tau) - op.tau;
        if (u_tau_p)
          new_mat = (*u_tau_p)[new_bl](dtau) * new_mat; // (interpolation)
        else {
          auto _ = triqs::arrays::range();
          for (int j = 0; j < dim; j++) new_mat(_, j) *= std::exp(-dtau * ad.get_eigenvalue(new_bl, j)); // Time-evolution
        }
        //std::cout << "new_mat 3: " << new_bl << " \n" << new_mat << "\n\n\n\n";
      }
      u_frame[new_bl] = new_mat;
    }
    return u_frame;
  }

} // namespace inchworm
