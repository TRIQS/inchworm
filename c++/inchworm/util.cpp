#include "./util.hpp"

#define NUM_BITS 10

namespace inchworm {

  uint64_t get_MSB(uint64_t a, int shift) {
    //std::printf("\nMSB %d ", shift);
    //print_binary(a, NUM_BITS);
    //print_binary(a >> shift, NUM_BITS - shift);
    return (a >> shift);
  }

  uint64_t get_LSB(uint64_t a, int shift) {
    //std::printf("\nLSB %d ", shift);
    //print_binary(a, NUM_BITS);
    //print_binary((a % (1 << shift)), shift);
    return (a % (1 << shift));
  }

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

  //
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

  //------------------------------

  many_body_operator create_effective_hyb(gf_struct_t const &gf_struct) {
    many_body_operator hyb_effective;
    for (auto const &[blname, blsize] : gf_struct) {
      for (auto [a, b] : product_range(blsize, blsize)) { hyb_effective += c_dag(blname, a) * c(blname, b); }
    }
    return hyb_effective;
  }

  std::pair<int, int> find_index(int number, std::vector<std::vector<fock_state_t>> fs) {
    for (int s = 0; s < fs.size(); s++) {
      for (int i = 0; i < fs[s].size(); i++) {
        if (number == fs[s][i]) return std::pair<int, int>(s, i);
      }
    }
    std::printf("error: number not found in the fock states\n");
    exit(1);
  }

  scalar_t trace(atom_diag const &ad_tot, std::function<double(double)> fct) {
    scalar_t trace_value = 0.0;
    auto es_full         = ad_tot.get_eigensystems();
    auto fs_full         = ad_tot.get_fock_states();

    for (int s = 0; s < ad_tot.n_subspaces(); s++)
      for (int i = 0; i < ad_tot.get_subspace_dim(s); i++) {
        scalar_t val = fct(es_full[s].eigenvalues[i] + ad_tot.get_gs_energy());
        trace_value += val;
        //std::printf("s=%d, i=%d, val= %f\n", s, i, val);
      }
    return trace_value;
  }

  frame_t partial_trace_bath(atom_diag const &ad_tot, atom_diag const &ad_imp, atom_diag const &ad_bath, double beta, double tau) {
    //TODO: incorporate in atom_diag and make it a member function: not possible anymore.

    for (int i = 0; i < (int)ad_imp.get_fops().data().size(); i++) {
      //std::printf("%d %d \n", i, (int)ad_imp.get_fops().data().size());
      if (ad_tot.get_fops().data()[i] != ad_imp.get_fops().data()[i]) {
        std::printf("error: the first indices of ad_tot should be the same as the one in ad_imp.\n");
        exit(1);
      }
    }

    //for (int i = 0; i < (int)ad_tot.get_fops().data().size(); i++) {
    //  std::cout << " " << ad_tot.get_fops().data()[i] << "\n";
    //}

    int linear_index = ad_imp.get_fops().data().size();
    //std::printf("li=%d\n",linear_index);

    frame_t u_frame_result = make_zero_propagator_frame(ad_imp);
    auto es_full           = ad_tot.get_eigensystems();
    auto fs_full           = ad_tot.get_fock_states();
    auto fs_loc            = ad_imp.get_fock_states();
    auto es_bath           = ad_bath.get_eigensystems();
    auto fs_bath           = ad_bath.get_fock_states();

    // printing:

    /*
    //print_eigensystems(ad_tot);
    for (int s = 0; s < ad_tot.n_subspaces(); s++) {
      for (int i = 0; i < ad_tot.get_subspace_dim(s); i++) {
        printf("  %2lu: ", fs_full[s][i]);
        print_binary(fs_full[s][i], ad_tot.get_fops().data().size());
      }
      std::cout << "\n";
    }
    std::cout << "\n";
    //print_eigensystems(ad_imp);
    for (int s = 0; s < ad_imp.n_subspaces(); s++) {
      for (int i = 0; i < ad_imp.get_subspace_dim(s); i++) {
        printf("  %2lu: ", fs_loc[s][i]);
        print_binary(fs_loc[s][i], ad_imp.get_fops().data().size());
      }
      std::cout << "\n";
    }
    //return 0.0;
    */

    for (int s = 0; s < ad_tot.n_subspaces(); s++) {
      EXPECTS(es_full[s].eigenvalues.size() == fs_full[s].size());
      int size      = ad_tot.get_subspace_dim(s);
      auto e_E_Udag = matrix_t{dagger(es_full[s].unitary_matrix)};

      for (int i = 0; i < size; i++) {
        //std::printf("-----> %lu   % 4.8f\n ", fs_full[s][i],  fct(es_full[s].eigenvalues[i] + ad_tot.get_gs_energy()));
        for (int j = 0; j < size; j++) e_E_Udag(i, j) *= std::exp(-tau * (es_full[s].eigenvalues[i] + ad_tot.get_gs_energy()));
      }
      //std::printf("\n");

      // e_H = exp[-tau * H] in the basis of the sites:
      auto e_H = es_full[s].unitary_matrix * e_E_Udag;

      /*
      // printing:
      std::printf("\nH[%d] %d \n   ", s);
      auto fss = ad_tot.get_fock_states()[s];
      for (int i = 0; i < ad_tot.get_subspace_dim(s); i++) {
        printf("  ");
        print_binary(fss[i], ad_tot.get_fops().data().size());
      }
    
      print_matrix(H);
      */

      for (int i = 0; i < size; i++) {
        uint64_t traced_idx1    = get_MSB(fs_full[s][i], linear_index); // the traced indices are the bath indices
        uint64_t preserved_idx1 = get_LSB(fs_full[s][i], linear_index); // the preserved indices are the impurity indices
        //std::printf("%u %u\n", traced_idx1, preserved_idx1);

        for (int j = 0; j < size; j++) {
          uint64_t traced_idx2    = get_MSB(fs_full[s][j], linear_index);
          uint64_t preserved_idx2 = get_LSB(fs_full[s][j], linear_index);
          //std::printf(" %u %u\n", traced_idx2, preserved_idx2);
          if (traced_idx1 == traced_idx2) {

            auto [s1, i1] = find_index(preserved_idx1, fs_loc);
            auto [s2, i2] = find_index(preserved_idx2, fs_loc);

            auto [s3, i3] = find_index(traced_idx1, fs_bath);
            //auto [s4, i4] = find_index(traced_idx2, fs_bath);

            if (s1 != s2) {
              std::printf("error: both block should be the same here\n");
              exit(1);
            }

            // the exponential function factor corresponds to e^[-(beta - tau) * H_bath]
            u_frame_result[s2](i1, i2) += e_H(i, j) * std::exp(-(beta - tau) * (es_bath[s3].eigenvalues[i3] + ad_bath.get_gs_energy()));
          }
        }
      }
    }

    // basis transformation to the atom_diag of the impurity:
    for (int s1 = 0; s1 < ad_imp.n_subspaces(); s1++) {
      /*
      std::printf("\nu[%d] \n   ", s1);
      auto fss1 = ad_imp.get_fock_states()[s1];
      for (int i = 0; i < ad_imp.get_subspace_dim(s1); i++) {
        printf("  ");
        print_binary(fss1[i], ad_imp.get_fops().data().size());
      }
      std::printf("\n");
      print_matrix(u_frame_result[s1]);
      */
      u_frame_result[s1] =
         dagger((ad_imp.get_eigensystems())[s1].unitary_matrix) * u_frame_result[s1] * ((ad_imp.get_eigensystems())[s1].unitary_matrix);
    }

    return u_frame_result;
  }
} // namespace inchworm
