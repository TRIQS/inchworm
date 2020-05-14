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

  void print_matrix(triqs::arrays::matrix<double> m) {

    for (int i = 0; i < first_dim(m); i++) {
      std::printf("\n [");
      for (int j = 0; j < second_dim(m); j++) {
        if (std::abs(m(i, j)) < 1e-12)
          std::printf("      .      ");
        else
          std::printf("% 13.6f", m(i, j));
      }
      std::printf("]");
    }
    std::printf("\n");

    //for (int i = 0; i < first_dim(m); i++) {
    //  for (int j = 0; j < second_dim(m); j++) { std::printf("% 5.6f ", m(i, j)); }
    //  std::printf("\n");
    //}
    //std::printf("\n\n");
  }

  void print_binary(unsigned int n, int total_bits) {
    if (n >= 0) {

      for (int i = 0; i < total_bits; i++) {
        int bit = (1 << (total_bits - i - 1));
        if ((n & bit) != 0)
          std::printf("1");
        else
          std::printf("0");
      }
      std::printf(" ");

    } else {
      std::printf(" negative binary? \n");
    }
  }

  /// Partial sum, tracing over indices above linear_index. Only the linear_index first degrees of freedom will be preserved.
  /**
     * @param ad atom_diag of the system considered here.
     * @param linear_index The linear index (i.e. number) of fundamental operator to be perserved, as defined by the fundamental operator set.
     * @param fct Function to be applied to eigenvalues in atom_diag.
     * @return The partial sum matrix of a function the Hamiltonian.
     */
  triqs::arrays::matrix<double> partial_sum(atom_diag const &ad, int linear_index, std::function<double(double)> fct) {
    //TODO: incorporate in atom_diag and make it a member function.
    int dim_partial = (1 << linear_index);
    //int dim_full    = ad.get_full_hilbert_space_dim();
    //int factor      = 1;//dim_full / dim_partial;
    //EXPECTS(dim_partial < dim_full);
    //EXPECTS(dim_full % dim_partial == 0);

    triqs::arrays::matrix<double> partial_sum(dim_partial, dim_partial);
    partial_sum = 0;

    auto es = ad.get_eigensystems();
    auto fs = ad.get_fock_states();
    EXPECTS(es.size() == fs.size());

    for (int s = 0; s < ad.n_subspaces(); s++) {
      EXPECTS(es[s].eigenvalues.size() == fs[s].size());
      int size    = ad.get_subspace_dim(s);
      auto E_Udag = dagger(es[s].unitary_matrix);
      for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++) E_Udag(i, j) *= fct(es[s].eigenvalues[i] + ad.get_gs_energy());
      auto H = es[s].unitary_matrix * E_Udag;

      for (int i = 0; i < size; i++) {
        uint64_t traced_idx1    = get_MSB(fs[s][i], linear_index);
        uint64_t preserved_idx1 = get_LSB(fs[s][i], linear_index);

        for (int j = 0; j < size; j++) {
          uint64_t traced_idx2    = get_MSB(fs[s][j], linear_index);
          uint64_t preserved_idx2 = get_LSB(fs[s][j], linear_index);
          if (traced_idx1 == traced_idx2) { partial_sum(preserved_idx1, preserved_idx2) += H(i, j); } // / factor; }
        }
      }
    }
    //std::printf("factor = %d, dim_full= %d, dim_partial= %d\n", factor, dim_full, dim_partial);
    return partial_sum;
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

  scalar_t trace(atom_diag const &ad_full, std::function<double(double)> fct) {
    scalar_t trace_value = 0.0;
    auto es_full         = ad_full.get_eigensystems();
    auto fs_full         = ad_full.get_fock_states();

    for (int s = 0; s < ad_full.n_subspaces(); s++)
      for (int i = 0; i < ad_full.get_subspace_dim(s); i++) {
        scalar_t val = fct(es_full[s].eigenvalues[i] + ad_full.get_gs_energy());
        trace_value += val;
        //std::printf("s=%d, i=%d, val= %f\n", s, i, val);
      }
    return trace_value;
  }

  u_frame_t partial_trace_bath(atom_diag const &ad_full, atom_diag const &ad_loc, atom_diag const &ad_bath, double beta, double dtau) {
    //TODO: incorporate in atom_diag and make it a member function.

    for (int i = 0; i < (int)ad_loc.get_fops().data().size(); i++) {
      //std::printf("%d %d \n", i, (int)ad_loc.get_fops().data().size());
      if (ad_full.get_fops().data()[i] != ad_loc.get_fops().data()[i]) {
        std::printf("error: the first indices of ad_full should be the same as the one in ad_loc.\n");
        exit(1);
      }

      //std::cout << ad_loc.get_fops().data()[i] << "\n";
    }

    //for (int i = 0; i < (int)ad_full.get_fops().data().size(); i++) {
    //  std::cout << " " << ad_full.get_fops().data()[i] << "\n";
    //}

    int linear_index = ad_loc.get_fops().data().size();
    //std::printf("li=%d\n",linear_index);

    u_frame_t u_frame_result = make_zero_propagator_frame(ad_loc);
    auto es_full             = ad_full.get_eigensystems();
    auto fs_full             = ad_full.get_fock_states();
    auto fs_target           = ad_loc.get_fock_states();
    auto es_bath             = ad_bath.get_eigensystems();
    auto fs_bath             = ad_bath.get_fock_states();

    //    auto fs_bath             = ad_bath.get_fock_states();
    /*
    //print_eigensystems(ad_full);
    for (int s = 0; s < ad_full.n_subspaces(); s++) {
      for (int i = 0; i < ad_full.get_subspace_dim(s); i++) {
        printf("  %2lu: ", fs_full[s][i]);
        print_binary(fs_full[s][i], ad_full.get_fops().data().size());
      }
      std::cout << "\n";
    }
    std::cout << "\n";
    //print_eigensystems(ad_loc);
    for (int s = 0; s < ad_loc.n_subspaces(); s++) {
      for (int i = 0; i < ad_loc.get_subspace_dim(s); i++) {
        printf("  %2lu: ", fs_target[s][i]);
        print_binary(fs_target[s][i], ad_loc.get_fops().data().size());
      }
      std::cout << "\n";
    }
    //return 0.0;
//*/

    for (int s = 0; s < ad_full.n_subspaces(); s++) {
      EXPECTS(es_full[s].eigenvalues.size() == fs_full[s].size());
      int size    = ad_full.get_subspace_dim(s);
      auto E_Udag = dagger(es_full[s].unitary_matrix);

      for (int i = 0; i < size; i++) {
        //std::printf("-----> %lu   % 4.8f\n ", fs_full[s][i],  fct(es_full[s].eigenvalues[i] + ad_full.get_gs_energy()));
        for (int j = 0; j < size; j++) E_Udag(i, j) *= std::exp(-dtau * (es_full[s].eigenvalues[i] + ad_full.get_gs_energy()));
        //fct(es_full[s].eigenvalues[i] + ad_full.get_gs_energy());
      }
      //std::printf("\n");
      auto H = es_full[s].unitary_matrix * E_Udag;
      //for (int i = 0; i < size; i++) {
      //  for (int j = 0; j < size; j++) H(i, j) *= std::exp(-dtau * (es_bath_full[s].eigenvalues[i] + ad_bath_full.get_gs_energy()));
      //}

      //print_matrix(H);
      for (int i = 0; i < size; i++) {
        uint64_t traced_idx1    = get_MSB(fs_full[s][i], linear_index); // the traced indices are the bath indices
        uint64_t preserved_idx1 = get_LSB(fs_full[s][i], linear_index); // the preserved indices are the impurity indices
        //std::printf("%u %u\n", traced_idx1, preserved_idx1);

        for (int j = 0; j < size; j++) {
          uint64_t traced_idx2    = get_MSB(fs_full[s][j], linear_index);
          uint64_t preserved_idx2 = get_LSB(fs_full[s][j], linear_index);
          //std::printf(" %u %u\n", traced_idx2, preserved_idx2);
          if (traced_idx1 == traced_idx2) {

            auto [s1, i1] = find_index(preserved_idx1, fs_target);
            auto [s2, i2] = find_index(preserved_idx2, fs_target);

            auto [s3, i3] = find_index(traced_idx1, fs_bath);
            auto [s4, i4] = find_index(traced_idx2, fs_bath);

            if (s1 != s2) {
              std::printf("error: both block should be the same here\n");
              exit(1);
            }

            // the exponential function factor corresponds to e^[-(beta - tau) * H_bath]
            u_frame_result[s1](i1, i2) += H(i, j) * std::exp(-(beta - dtau) * (es_bath[s3].eigenvalues[i3] + ad_bath.get_gs_energy()));

          } // partial_sum(preserved_idx1, preserved_idx2) += H(i, j); }
        }
      }
    }

    // basis transformation to the atom_diag of the impurity
    for (int s1 = 0; s1 < ad_loc.n_subspaces(); s1++) {
      u_frame_result[s1] =
         ((ad_loc.get_eigensystems())[s1].unitary_matrix * u_frame_result[s1]) * dagger((ad_loc.get_eigensystems())[s1].unitary_matrix);
    }

    return u_frame_result;
  }
} // namespace inchworm
