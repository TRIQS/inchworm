#include "./util.hpp"

namespace inchworm {

  void print_energies(std::vector<std::vector<double>> const &E) {
    for (auto sp : E) {
      for (auto l : sp) { std::printf("% 2.3f ", l); }
      std::printf("\n");
    }
    std::printf("\n");
  }

  void print_eigensystems(atom_diag const &ad) {
    for (auto sp : ad.get_eigensystems()) {
      for (auto l : sp.eigenvalues) { std::printf("% 2.3f ", l); }
      std::printf("\n\n");

      for (int i = 0; i < sp.eigenvalues.size(); i++) {
        for (int j = 0; j < sp.eigenvalues.size(); j++) { std::printf("% 2.3f ", sp.unitary_matrix(i, j)); }
        std::printf("\n");
      }
      //for (auto u : sp.unitary_matrix) { TRIQS_PRINT(u); }
      std::printf("\n\n");
    }
    std::printf("\n");
  }

  void print_binary(unsigned int n, int total_bits) {

    for (int i = 0; i < total_bits; i++) {
      int bit = (1 << (total_bits - i - 1));
      if ((n & bit) != 0)
        std::printf("1");
      else
        std::printf("0");
    }
    std::printf(" ");
  }

  void print_atom_diag(atom_diag const &ad) {

    int N = 0;
    for (int s = 0; s < ad.n_subspaces(); s++) {
      auto sp     = ad.get_eigensystems()[s];
      auto E_Udag = matrix_t{dagger(sp.unitary_matrix)};
      for (int i = 0; i < sp.eigenvalues.size(); i++) {
        for (int j = 0; j < sp.eigenvalues.size(); j++) E_Udag(i, j) *= sp.eigenvalues[i] + ad.get_gs_energy();
      }
      auto H = sp.unitary_matrix * E_Udag;

      std::printf("\nblock: %d \n   ", N++);
      auto fss = ad.get_fock_states()[s];
      for (int i = 0; i < ad.get_subspace_dim(s); i++) {
        printf("  ");
        print_binary(fss[i], ad.get_fops().data().size());
      }
      std::printf("\n");
      print_matrix(H);
    }
    std::printf("\n");
  }

  void print_matrix(triqs::arrays::matrix<double> m, double factor) {

    for (int i = 0; i < m.shape()[0]; i++) {
      std::printf("\n [");
      for (int j = 0; j < m.shape()[1]; j++) {
        if (std::abs(m(i, j)) == 0.0)
          std::printf("  .        ");
        else
          std::printf("% 10.3e ", m(i, j) * factor);
      }
      std::printf("]");
    }
    std::printf("\n");
  }

  //
  void fprint(u_tau_t u_tau, int N_tau) {
    FILE *f = fopen("u_tau.dat", "w");
    for (int i_tau = 0; i_tau < N_tau; i_tau++) {
      fprintf(f, "%d  ", i_tau);
      for (int bl = 0; bl < u_tau.size(); bl++) print_matrix((matrix_t)u_tau[bl][i_tau]);
      fprintf(f, "\n");
    }
    fclose(f);
    return;
  }

  //
  void print(u_partial_t u_partial) {
    for (auto &[bl, mat] : u_partial) {
      std::printf("\n\nblock: %d\n", bl);
      print_matrix(mat);
    }
    std::cout << "\n";
    return;
  }

  //
  void print(frame_t u_frame) {
    for (int bl = 0; bl < u_frame.size(); bl++) {
      std::printf("\n\nblock: %d\n", bl);
      print_matrix(u_frame[bl]);
    }
    std::cout << "\n";
    return;
  }

  //
  void print(u_tau_t u_tau, double tau) {
    for (int bl = 0; bl < u_tau.size(); bl++) { print_matrix((matrix_t)u_tau[bl](tau)); }
    std::cout << "\n";
    return;
  }

  //
  void print(u_tau_t u_tau, int frame_number) {
    for (int bl = 0; bl < u_tau.size(); bl++) { std::cout << std::setprecision(10) << u_tau[bl][frame_number]; }
    std::cout << "\n";
    return;
  }

  //------------------------------

  many_body_operator create_effective_hyb(gf_struct_t const &gf_struct) {
    many_body_operator hyb_effective;
    for (auto const &[blname, blsize] : gf_struct) {
      for (auto [a, b] : product_range(blsize, blsize)) { hyb_effective += c_dag(blname, a) * c(blname, b); }
    }
    return hyb_effective;
  }

  //------------------------------

  /**
   * Split the Fockstate (bitset) at a given bit index
   *
   * Example:
   * (1 1 1 0 1 0 0..) -> [(1 1 1 0 0..), (0 1 0 0...)]
   *       ^ bit_index = 3
   */
  std::pair<uint64_t, uint64_t> split_fs(uint64_t fs, int bit_index) {
    //std::cout << std::bitset<10>(fs) << "\n";
    uint64_t fs_left  = fs % (1 << bit_index); // All bits left of bit_index
    uint64_t fs_right = fs >> bit_index;       // All bits right of bit_index
    return {fs_left, fs_right};
  }

  /**
   * Given a Fockstate, return the block and index
   * for the associated atom_diag object
   */
  std::pair<int, int> fs_to_bl_and_idx(uint64_t fs, atom_diag const &ad) {
    auto const &all_fs = ad.get_fock_states();
    for (int bl = 0; bl < all_fs.size(); bl++) {
      for (int idx = 0; idx < all_fs[bl].size(); idx++) {
        if (fs == all_fs[bl][idx]) return {bl, idx};
      }
    }
    std::printf("error: number not found in the fock states\n");
    std::abort();
  }

  frame_t partial_trace_bath(atom_diag const &ad_tot, atom_diag const &ad_imp, atom_diag const &ad_bath, double beta, double tau) {

    for (auto i : range(ad_imp.get_fops().size())) {
      if (ad_tot.get_fops().data()[i] != ad_imp.get_fops().data()[i]) {
        std::printf("error: the first indices of ad_tot should be the same as the one in ad_imp.\n");
        std::abort();
      }
    }

    // Given atom_diag object, calculate e^[-tau * H] in Fockstate Basis
    auto calc_e_H_fs = [](atom_diag const &ad, double tau) {
      // e^[-tau * H] in eigenbasis
      auto e_H_tau = make_bare_u_frame(ad, tau);

      // Rotate to Fockstate Basis
      auto e_H_tau_fs = e_H_tau;
      for (auto bl : range(ad.n_subspaces())) {
        auto rot       = ad.get_eigensystems()[bl].unitary_matrix;
        e_H_tau_fs[bl] = rot * e_H_tau[bl] * dagger(rot);
      }
      return e_H_tau_fs;
    };
    auto e_H_tot_tau_fs             = calc_e_H_fs(ad_tot, tau);
    auto e_H_bath_beta_minus_tau_fs = calc_e_H_fs(ad_bath, beta - tau);

    // The result container
    frame_t utau_imp = make_frame(ad_imp.get_subspace_dims());

    auto const &all_fs_tot = ad_tot.get_fock_states();

    for (auto [bl, bl_size] : enumerate(ad_tot.get_subspace_dims())) {

      for (auto [i, j] : product_range(bl_size, bl_size)) {

        int nfops_imp              = ad_imp.get_fops().size();
        auto [fs_imp_i, fs_bath_i] = split_fs(all_fs_tot[bl][i], nfops_imp);
        auto [fs_imp_j, fs_bath_j] = split_fs(all_fs_tot[bl][j], nfops_imp);

        if (fs_bath_i == fs_bath_j) {

          auto [bl_imp, i_imp]   = fs_to_bl_and_idx(fs_imp_i, ad_imp);
          auto [bl_imp_j, j_imp] = fs_to_bl_and_idx(fs_imp_j, ad_imp);
          ASSERT(bl_imp == bl_imp_j);

          auto [bl_bath, i_bath] = fs_to_bl_and_idx(fs_bath_i, ad_bath);

          utau_imp[bl_imp](i_imp, j_imp) += e_H_tot_tau_fs[bl](i, j) * e_H_bath_beta_minus_tau_fs[bl_bath](i_bath, i_bath);
        }
      }
    }

    // Rotate to impurity eigenbasis and normalize by Z_bath
    auto Z_bath = trace(make_bare_u_frame(ad_bath, beta));
    for (auto bl_imp : range(ad_imp.n_subspaces())) {
      auto rot = ad_imp.get_eigensystems()[bl_imp].unitary_matrix;
      utau_imp[bl_imp] = 1.0 / Z_bath * dagger(rot) * utau_imp[bl_imp] * rot;
    }

    return utau_imp;
  }
} // namespace inchworm
