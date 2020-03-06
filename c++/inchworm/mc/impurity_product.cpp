#include "./impurity_product.hpp"

namespace inchworm {

  u_tau_t make_propagator(atom_diag const &h_diag, double beta, int n_tau) {
    // this create the propagator and assign identity to the first frame (or time) of the propagator.
    int n_sub = h_diag.n_subspaces();
    triqs::hilbert_space::gf_struct_t propagator_struct;

    for (int i = 0; i < n_sub; i++) {
      std::vector<std::variant<int, std::string>> l(h_diag.get_subspace_dim(i));
      std::iota(l.begin(), l.end(), 0);
      propagator_struct.push_back(std::make_pair(std::to_string(i), l));
    }

    auto u_tau = u_tau_t{{beta, Fermion, n_tau}, propagator_struct};
    for (auto &block : u_tau) block[0] = 1; // identity at time zero
    return u_tau;
  }

  //
  void print(u_tau_t u_tau, double tau) {
    for (int bl = 0; bl < u_tau.size(); bl++) { std::cout << u_tau[bl](tau); }
    std::cout << "\n";
    return;
  }

  //
  void print(u_tau_t u_tau, int frame_number) {
    for (int bl = 0; bl < u_tau.size(); bl++) { std::cout << u_tau[bl][frame_number]; }
    std::cout << "\n";
    return;
  }

  //
  void assign_u_frame_to_propagator(u_tau_t &u_tau, u_frame_t const &u_frame, int frame_number, scalar_t factor) {
    for (int bl = 0; bl < u_tau.size(); bl++) u_tau[bl][frame_number] = factor * u_frame[bl];
    return;
  }

  u_tau_t make_ED_propagator(atom_diag const &ad_tot, atom_diag const &ad_atom, atom_diag const &ad_bath, double beta, int n_tau) {

    auto E0       = ad_tot.get_gs_energy();
    u_tau_t u_tau = make_propagator(ad_atom, beta, n_tau);

    for (int i_tau = 0; i_tau < n_tau; i_tau++) {
      double dtau  = beta * i_tau / (n_tau - 1);
      auto u_frame = partial_trace(ad_tot, ad_atom, [dtau, E0](double E) { return std::exp(-dtau * (E - E0)); });
      auto Z_bath  = trace(ad_bath, [dtau, E0](double E) { return std::exp(-dtau * (E - E0)); });
      assign_u_frame_to_propagator(u_tau, u_frame, i_tau, 1. / Z_bath);
    }

    return u_tau;
  }

  //
  u_frame_t propagator_product(atom_diag const &ad, time_diagram_t const &diagram, double tau, u_tau_t const *const u_tau_p) {

    if (diagram.size() == 0) return make_bare_propagator_frame(ad, tau, false);

    constexpr bool set_gs_to_0 = false;
    u_frame_t u_frame          = make_zero_propagator_frame(ad);
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
          new_mat(j, j) = std::exp(
             -dtau * (ad.get_eigenvalue(initial_bl, j) + (set_gs_to_0 ? 0. : ad.get_gs_energy()))); // Create time-evolution matrix e^-H(tau-tau_max)
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
          for (int j = 0; j < ad.get_subspace_dim(new_bl); j++) {
            //std::printf("new_bl %d, j %d \n", new_bl, j);
            //std::printf("ad.get_subspace_dim(new_bl) = %d \n", ad.get_subspace_dim(new_bl));
            new_mat(j, _) *=                                                                               // ATTENTION!
               std::exp(-dtau * (ad.get_eigenvalue(new_bl, j) + (set_gs_to_0 ? 0. : ad.get_gs_energy()))); // bare imaginary time evolution
          }
        }
        //std::cout << "new_mat 3: " << new_bl << " \n" << new_mat << "\n\n\n\n";
      }
      u_frame[new_bl] = new_mat;
    }
    return u_frame;
  }

} // namespace inchworm
