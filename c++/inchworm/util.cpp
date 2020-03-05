#include "./util.hpp"

namespace inchworm {

  /// Partial sum, tracing over indices above linear_index. Only the linear_index first degrees of freedom will be preserved.
  /**
     * @param ad atom_diag of the system considered here.
     * @param linear_index The linear index (i.e. number) of fundamental operator to be perserved, as defined by the fundamental operator set.
     * @param fct Function to be applied to eigenvalues in atom_diag.
     * @return The partial sum matrix of a function the Hamiltonian.
     */
  triqs::arrays::matrix<double> partial_sum(triqs::atom_diag::atom_diag<false> const &ad, int linear_index, std::function<double(double)> fct) {
    //TODO: incorporate in atom_diag and make it a member function.
    int dim_partial = (1 << linear_index);
    int dim_full    = ad.get_full_hilbert_space_dim();
    //int factor      = 1;//dim_full / dim_partial;
    EXPECTS(dim_partial < dim_full);
    EXPECTS(dim_full % dim_partial == 0);

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

  triqs::arrays::matrix<double> partial_sum2(triqs::atom_diag::atom_diag<false> const &ad, int linear_index, std::function<double(double)> fct) {
    //TODO: incorporate in atom_diag and make it a member function.
    int dim_partial = (1 << linear_index);
    int dim_full    = ad.get_full_hilbert_space_dim();
    //int factor      = 1;//dim_full / dim_partial;
    EXPECTS(dim_partial < dim_full);
    EXPECTS(dim_full % dim_partial == 0);

    triqs::arrays::matrix<double> partial_sum(dim_partial, dim_partial);
    partial_sum = 0;

    auto es = ad.get_eigensystems();
    auto fs = ad.get_fock_states();
    EXPECTS(es.size() == fs.size());

    for (int s = 0; s < ad.n_subspaces(); s++) {
      EXPECTS(es[s].eigenvalues.size() == fs[s].size());
      int size    = ad.get_subspace_dim(s);
      auto E_Udag = dagger(es[s].unitary_matrix);
      for (int i = 0; i < size; i++) {
        std::printf("%lu ", fs[s][i]);
        for (int j = 0; j < size; j++) E_Udag(i, j) *= fct(es[s].eigenvalues[i] + ad.get_gs_energy());
      }
      std::printf("\n");
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

} // namespace inchworm
