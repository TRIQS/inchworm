#include "./impurity_product.hpp"

namespace inchworm {

  u_partial_t impurity_product(atom_diag const &ad, time_diagram_t const &diagram, double tau_max, double tau_min, u_tau_t const *const u_tau_p) {

    EXPECTS(tau_min < tau_max);

    std::vector<int> op_idx;
    for (int i = 0; (i < diagram.size()); i++) {
      auto const &op = diagram.op_list[i];
      if ((op.tau >= tau_min) and (op.tau <= tau_max)) op_idx.push_back(i);
    }

    // Treat the trivial case of zero operators separately
    if (op_idx.size() == 0) {
      if (u_tau_p) {
        frame_t u_frame = make_zero_frame(ad.get_subspace_dims());
        for (int bl = 0; bl < ad.n_subspaces(); bl++) {
          u_frame[bl] = (*u_tau_p)[bl](tau_max - tau_min); // (interpolation)
        }
        return make_u_partial(u_frame);
      } else {
        return make_u_partial(make_bare_u_frame(ad, tau_max - tau_min));
      }
    }

    u_partial_t u_partial(ad.n_subspaces());

    // loop over all the block of the propagator. initial_bl is the starting block before any d/d_dag operator.
    // new_bl is the bloc after apply "i" operator d/d_dag. After i=diagram.size() appication of operator (d/d_dag)
    // new_bl is the last block onto which initial_bl is mapped. We first determine what is this final block and
    // if it is not -1, we proceed to calculate matrix multiplications.
    for (int initial_bl = 0; initial_bl < ad.n_subspaces(); initial_bl++) {
      int dim      = ad.get_subspace_dim(initial_bl);
      int new_bl   = initial_bl;
      auto new_mat = matrix_t{};

      // first calculate the final block withou matrix multiplication:
      for (auto i : op_idx) {
        //std::printf("test %d,  %d\n", i, diagram.size());
        if (new_bl == -1) break;
        auto const &op = diagram.op_list[i];
        // apply d or d_dag operator. Note that the d/d_dag operator in our formalism corresponds to c/cdag operator of atom_diag
        new_bl = (op.dag ? ad.cdag_connection(op.linear_index, new_bl) : ad.c_connection(op.linear_index, new_bl));
      }
      //std::printf("\n");
      //std::printf("bloc: %d, goes to: %d \n", initial_bl, new_bl);
      if (new_bl == -1) {
        u_partial[initial_bl] = {-1, {}};
        continue; // do not proceed to matrix multipication if the final bloc is -1.
      }
      new_bl = initial_bl;

      // start by calculating U(tau_0)
      double dtau = diagram.op_list[op_idx[0]].tau - tau_min;
      //double dtau2 = 0.0;
      if (u_tau_p) {                            // if u_tau is defined, calculate inchworm case
        new_mat = (*u_tau_p)[initial_bl](dtau); // (interpolation)
      } else {                                  // if u_tau is not defined, calculate cthyb case
        new_mat = matrix_t(dim, dim);           //zeros?
        new_mat = 0;
        for (int j = 0; j < dim; j++)
          new_mat(j, j) = std::exp(-dtau * (ad.get_eigenvalue(initial_bl, j) + (ad.get_gs_energy()))); // Create time-evolution matrix e^(-H*tau)
      }

      // matrix multiplication = U(tau-tau_2k)*D_2k*U(tau_2k-tau_2k-1)*...*U(tau_2-tau_1)*D_1*U(tau_1-tau_0)*D_0*U(tau_0)
      // at the beginning of the loop, new_mat = U(tau_0)
      for (auto i : op_idx) {
        auto const &op = diagram.op_list[i];

        // apply d or d_dag operator. Note that the d/d_dag operator in our formalism corresponds to c/cdag operator of atom_diag
        new_mat = (op.dag ? ad.cdag_matrix(op.linear_index, new_bl) : ad.c_matrix(op.linear_index, new_bl)) * new_mat;
        new_bl  = (op.dag ? ad.cdag_connection(op.linear_index, new_bl) : ad.c_connection(op.linear_index, new_bl));

        if (u_tau_p) {

          //dtau2 = 0.0;
          if (i == op_idx.back()) { // if last point of the diagram
            dtau = tau_max - op.tau;
          } else {                                      // if not last point of the diagram
            dtau = diagram.op_list[i + 1].tau - op.tau; //FIXME: i+1: i do not like it.
          }

          new_mat = (*u_tau_p)[new_bl](dtau) * new_mat; // (interpolation)
        } else {                                        // use bare propagator (cthyb)
          auto _ = triqs::arrays::range();
          dtau   = (i == (op_idx.back()) ? tau_max : diagram.op_list[i + 1].tau) - op.tau;

          for (int j = 0; j < ad.get_subspace_dim(new_bl); j++) {
            new_mat(j, _) *= std::exp(-dtau * (ad.get_eigenvalue(new_bl, j) + (ad.get_gs_energy()))); // bare imaginary time evolution
          }
        }
      }
      u_partial[initial_bl] = {new_bl, new_mat};
    }
    return u_partial;
  }

} // namespace inchworm
