#include "./u_frame.hpp"

namespace inchworm {

  frame_t make_u_frame(u_partial_t const &up) {
    frame_t res;
    for (int i = 0; i < up.size(); ++i) {
      auto &[bl, mat] = up[i];
      EXPECTS(i == bl || bl == -1);
      res.emplace_back(mat);
    }
    return res;
  }

  u_partial_t make_u_partial(frame_t const &u) {
    u_partial_t res;
    for (int i = 0; i < u.size(); ++i) { res.emplace_back(i, u[i]); }
    return res;
  }

  frame_t make_zero_propagator_frame(atom_diag const &ad) {
    frame_t u_frame(ad.n_subspaces());

    for (int bl = 0; bl < ad.n_subspaces(); bl++) {
      u_frame[bl] = matrix_t(ad.get_subspace_dim(bl), ad.get_subspace_dim(bl)); //  use zeros<> ?? check
      u_frame[bl] = 0;
    }
    return u_frame;
  }

  frame_t make_bare_propagator_frame(atom_diag const &ad, double tau, bool set_gs_to_0) {
    frame_t u_frame(ad.n_subspaces());

    for (int bl = 0; bl < ad.n_subspaces(); bl++) {
      int dim     = ad.get_subspace_dim(bl);
      u_frame[bl] = matrix_t(dim, dim); //  use zeros<> ?? check
      u_frame[bl] = 0;
      for (int j = 0; j < dim; j++) u_frame[bl](j, j) = std::exp(-tau * (ad.get_eigenvalue(bl, j) + (set_gs_to_0 ? 0. : ad.get_gs_energy())));
    }
    return u_frame;
  }

  frame_t make_bare_g_frame(atom_diag const &ad_imp, u_tau_t const &u_tau, std::map<int, std::pair<int, int>> const &map_lin_idx_to_block_inner,
                            gf_struct_t const &gf_struct, double tau_split, double beta) {
    u_partial_t l, r;

    for (int i = 0; i < u_tau.size(); ++i) {
      l.push_back({i, u_tau[i](beta - tau_split)});
      r.push_back({i, u_tau[i](tau_split)});
    }

    return make_g_frame_from_l_and_r(ad_imp, map_lin_idx_to_block_inner, gf_struct, l, r);
  }

  frame_t make_frame(std::vector<long> const &shape_of_frame) {
    auto res = frame_t{};

    for (auto n : shape_of_frame) {
      res.emplace_back(matrix_t(n, n));
      res.back() = 0.;
    }
    return res;
  }

  // Create an empty frame (block diagonal matrix: vector of matrix_t)
  frame_t make_frame(gf_struct_t const &gf_struct) {
    auto res = frame_t{};

    for (auto bl : range(gf_struct.size())) {
      auto &[blname, blsize] = gf_struct[bl];
      res[bl] = matrix_t{blsize, blsize};
      res[bl] = 0.;
    }
    return res;
  }

  u_partial_t operator*(u_partial_t const &l, u_partial_t const &r) {

    EXPECTS(l.size() == r.size());

    auto res = u_partial_t{};
    for (int i = 0; i < l.size(); ++i) {
      // u[bl] = up[bl, blp] * up[blp, bl]
      // l_bl <- r_bl <- i
      auto &[r_bl, r_mat] = r[i];
      if (r_bl == -1)
        res.emplace_back(-1, matrix_t{});
      else {
        auto &[l_bl, l_mat] = l[r_bl];
        if (l_bl == -1)
          res.emplace_back(-1, matrix_t{});
        else
          res.emplace_back(l_bl, l_mat * r_mat);
      }
    }
    return res;
  }

  u_partial_t apply_op_from_right(u_partial_t const &l, int lin_index, bool op_dag, atom_diag const &ad) {

    auto res = u_partial_t{};
    for (int i = 0; i < l.size(); ++i) {
      // u[bl] = up[bl, blp] * d_dag
      // l_bl <- r_bl <- i
      //auto &[r_bl, r_mat] = r[i];

      auto r_bl  = (op_dag ? ad.cdag_connection(lin_index, i) : ad.c_connection(lin_index, i));
      auto r_mat = (op_dag ? ad.cdag_matrix(lin_index, i) : ad.c_matrix(lin_index, i));

      if (r_bl == -1) {
        res.emplace_back(-1, matrix_t{});
      } else {
        auto &[l_bl, l_mat] = l[r_bl];
        if (l_bl == -1) {
          res.emplace_back(-1, matrix_t{});
        } else {
          res.emplace_back(l_bl, l_mat * r_mat);
        }
      }
    }
    return res;

  } // namespace inchworm

  double frobenius_norm(frame_t const &g_frame) {
    double val = 0;
    for (auto const &mat : g_frame) {
      double norm = frobenius_norm(mat);
      val += norm * norm;
      //double norm = trace(B);
      //val += norm ;
    }
    return std::sqrt(val);
    //return val;
  }

  double frobenius_norm(u_partial_t const &u_partial) {
    double val = 0;
    for (auto const &[bl, mat] : u_partial) {
      double norm = frobenius_norm(mat);
      val += norm * norm;
      //double norm = trace(B);
      //val += norm ;
    }
    return std::sqrt(val);
    //return val;
  }

  // calculate the trace of the u_frame block diagonal matrix:
  double trace(frame_t const &u_frame) {
    double val = 0;
    for (auto const &B : u_frame) val += trace(B);
    return val;
  }

  double relative_distance(frame_t const &l, frame_t const &r) {

    frame_t diff = l;
    for (auto bl : range(l.size())) { diff[bl] = l[bl] - r[bl]; }

    auto norm_l         = frobenius_norm(l);
    auto norm_r         = frobenius_norm(r);
    auto norm_l_minus_r = frobenius_norm(diff);

    return norm_l_minus_r / std::max(norm_l, norm_r);
  }

  frame_t get_frame(u_tau_t const & u_tau, int idx){
    frame_t res;
    for(auto ubl: u_tau)
      res.push_back(ubl[idx]);
    return res;
  }

  double relative_distance(u_tau_t const &l, u_tau_t const &r) {
    double dist = 0.0;
    for(int i = 0; i < l[0].mesh().size(); ++i)
      dist = std::max(dist, relative_distance(get_frame(l,i), get_frame(r,i))); 
    return dist;
  }

  void print(frame_t const &u_frame, double factor) {
    for (auto block : u_frame) {
      block *= factor;
      std::cout << block;
    }
    std::cout << "\n";
  }

  frame_t make_g_frame_from_l_and_r(atom_diag const &ad_imp, std::map<int, std::pair<int, int>> const &map_lin_idx_to_block_inner,
                                    gf_struct_t const &gf_struct, u_partial_t const &l, u_partial_t const &r) {

    frame_t g_frame = make_frame(gf_struct);

    // G[g_bl][tau][in,in_dag] = -<T c[g_bl][in](tau) cdag[g_bl][in_dag](0)>
    // i ~= g_bl     + in
    // j ~= g_bl_dag + in_dag
    int n_ops = ad_imp.get_fops().data().size();
    for (int i = 0; i < n_ops; ++i) {
      auto [g_bl, in] = map_lin_idx_to_block_inner.at(i);
      for (int j = 0; j < n_ops; ++j) {
        auto [g_bl_dag, in_dag] = map_lin_idx_to_block_inner.at(j);
        if (g_bl != g_bl_dag) continue;

        // r * ddag_j[bl_idx1](0)
        u_partial_t rddag = apply_op_from_right(r, j, true, ad_imp);

        // l * d_i[bl_idx2](tau)
        u_partial_t ld = apply_op_from_right(l, i, false, ad_imp);

        auto prod = make_u_frame(ld * rddag);

        for (int bl0 = 0; bl0 < ad_imp.n_subspaces(); ++bl0) {
          g_frame[g_bl](in, in_dag) -= trace(prod[bl0]); //FIXME check the order of in and in_dag to be sure.
        }
      }
    }

    return g_frame;
  }

} // namespace inchworm
