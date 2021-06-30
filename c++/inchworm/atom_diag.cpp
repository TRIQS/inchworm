#include "./atom_diag.hpp"
#include "./u_frame.hpp"

namespace inchworm {

  many_body_operator create_effective_hyb(gf_struct_t const &gf_struct) {
    many_body_operator hyb_effective;
    for (auto const &[blname, blsize] : gf_struct) {
      for (auto [a, b] : product_range(blsize, blsize)) { hyb_effective += c_dag(blname, a) * c(blname, b); }
    }
    return hyb_effective;
  }

  frame_t make_g_frame_from_l_and_r(atom_diag const &ad_imp, gf_struct_t const &gf_struct, u_partial_t const &l, u_partial_t const &r) {

    frame_t g_frame = make_zero_frame(gf_struct);

    // G[bl][tau][i,j] = -<T c[bl][i](tau) cdag[bl][j](0)>
    for (int bl : range(gf_struct.size())) {
      auto const &[bl_name, bl_size] = gf_struct[bl];

      for (auto [i, j] : product_range(bl_size, bl_size)) {

        auto l_x_di    = l * get_op_block_matrix(ad_imp, bl_name, i, false);
        auto r_x_djdag = r * get_op_block_matrix(ad_imp, bl_name, j, true);
        auto prod      = make_frame(l_x_di * r_x_djdag);

        g_frame[bl](i, j) -= trace(prod);
      }
    }

    return g_frame;
  }

  frame_t make_bare_u_frame(atom_diag const &ad, double tau) {
    auto u_frame = make_zero_frame(ad.get_subspace_dims());
    for (auto [bl, bl_size] : enumerate(ad.get_subspace_dims()))
      for (int i : range(bl_size)) u_frame[bl](i, i) = std::exp(-tau * ad.get_eigenvalue(bl, i));
    return u_frame;
  }

  u_tau_t make_ED_propagator(atom_diag const &ad_tot, atom_diag const &ad_imp, atom_diag const &ad_bath, double beta, int n_tau) {

    auto u_tau = u_tau_t{{beta, Fermion, n_tau}, ad_imp.get_subspace_dims()};

    double dtau = beta / (n_tau - 1.);
    for (int i_tau = 0; i_tau < n_tau; i_tau++) {
      auto u_frame = partial_trace_bath(ad_tot, ad_imp, ad_bath, beta, dtau * i_tau);
      set_frame(u_frame, u_tau, i_tau);
    }

    return u_tau;
  }

  frame_t make_bare_g_frame(atom_diag const &ad_imp, u_tau_t const &u_tau, gf_struct_t const &gf_struct, double tau_split, double beta) {
    u_partial_t l(u_tau.size()), r(u_tau.size());

    for (int i = 0; i < u_tau.size(); ++i) {
      l[i] = {i, u_tau[i](beta - tau_split)};
      r[i] = {i, u_tau[i](tau_split)};
    }

    return make_g_frame_from_l_and_r(ad_imp, gf_struct, l, r);
  }

  u_partial_t get_op_block_matrix(atom_diag const &ad, std::string const &bl_name, int idx, bool op_dag) {
    auto res = u_partial_t(ad.n_subspaces());

    auto lidx = ad.get_fops()[{bl_name, idx}];
    for (auto bl_in : range(ad.n_subspaces())) {
      auto bl_out = (op_dag ? ad.cdag_connection(lidx, bl_in) : ad.c_connection(lidx, bl_in));
      auto matrix = (op_dag ? ad.cdag_matrix(lidx, bl_in) : ad.c_matrix(lidx, bl_in));
      res[bl_in]  = {bl_out, matrix};
    }

    return res;
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
    frame_t utau_imp = make_zero_frame(ad_imp.get_subspace_dims());

    auto const &all_fs_tot = ad_tot.get_fock_states();

    for (auto [bl, bl_size] : enumerate(ad_tot.get_subspace_dims())) {

      for (auto [i, j] : product_range(bl_size, bl_size)) {

        int nfops_imp              = ad_imp.get_fops().size();
        auto [fs_imp_i, fs_bath_i] = split_fs(all_fs_tot[bl][i], nfops_imp);
        auto [fs_imp_j, fs_bath_j] = split_fs(all_fs_tot[bl][j], nfops_imp);

        if (fs_bath_i == fs_bath_j) {

          auto [bl_imp_i, i_imp] = fs_to_bl_and_idx(fs_imp_i, ad_imp);
          auto [bl_imp_j, j_imp] = fs_to_bl_and_idx(fs_imp_j, ad_imp);

          if (bl_imp_i != bl_imp_j) TRIQS_RUNTIME_ERROR << "Block structore of ad_imp incompatible with hybridization";

          auto [bl_bath, i_bath] = fs_to_bl_and_idx(fs_bath_i, ad_bath);

          utau_imp[bl_imp_i](i_imp, j_imp) += e_H_tot_tau_fs[bl](i, j) * e_H_bath_beta_minus_tau_fs[bl_bath](i_bath, i_bath);
        }
      }
    }

    // Rotate to impurity eigenbasis and normalize by Z_bath
    auto Z_bath = trace(make_bare_u_frame(ad_bath, beta));
    for (auto bl_imp : range(ad_imp.n_subspaces())) {
      auto rot         = ad_imp.get_eigensystems()[bl_imp].unitary_matrix;
      utau_imp[bl_imp] = 1.0 / Z_bath * dagger(rot) * utau_imp[bl_imp] * rot;
    }

    // Correct for using different reference energies in imp, bath and tot
    double factor = std::exp((ad_bath.get_gs_energy() + ad_imp.get_gs_energy() - ad_tot.get_gs_energy()) * tau);

    return factor * utau_imp;
  }

} // namespace inchworm
