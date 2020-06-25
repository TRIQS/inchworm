#include "./u_frame.hpp"

namespace inchworm {
  // Create an empty frame (block diagonal matrix: vector of matrix_t)
  frame_t make_zero_propagator_frame(atom_diag const &ad) {
    frame_t u_frame(ad.n_subspaces());

    for (int bl = 0; bl < ad.n_subspaces(); bl++) {
      u_frame[bl] = matrix_t(ad.get_subspace_dim(bl), ad.get_subspace_dim(bl)); //  use zeros<> ?? check
      u_frame[bl] = 0;
    }
    return u_frame;
  }

  // initialize bare propagator frame U_0 = exp(-tau H_loc) in the diagonal basis of H_loc
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

  u_partial_t make_u_partial(frame_t const &u) {
    u_partial_t res;
    for (int i = 0; i < u.size(); ++i) { res.emplace_back(i, u[i]); }
    return res;
  }

  frame_t make_u_frame(u_partial_t const &up) {
    frame_t res;
    for (int i = 0; i < up.size(); ++i) {
      auto &[bl, mat] = up[i];
      EXPECTS(i == bl || bl == -1);
      res.emplace_back(mat);
    }
    return res;
  }

  frame_t make_frame(std::vector<long> const & shape_of_frame) {
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

    for (auto const &[bl, idxlst] : gf_struct) {
      res.emplace_back(idxlst.size(), idxlst.size());
      res.back() = 0.;
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

  //std::vector<std::pair<int /*op.linear_index*/, matrix_t /*c_matrix or cdag_matrix*/>> c_lst;
  //c_lst[bl][lidx]

  u_partial_t apply_op_from_right(u_partial_t const &l, int lin_index, bool op_dag, atom_diag const &ad) {
    auto res = u_partial_t{};
    for (int i = 0; i < l.size(); ++i) {
      // u[bl] = up[bl, blp] * d_dag
      // l_bl <- r_bl <- i
      //auto &[r_bl, r_mat] = r[i];

      auto r_bl = (op_dag ? ad.cdag_connection(lin_index, i) : ad.c_connection(lin_index, i));
      if (r_bl == -1) {
        res.emplace_back(-1, matrix_t{});
      } else {
        auto &[l_bl, l_mat] = l[r_bl];
        if (l_bl == -1) {
          res.emplace_back(-1, matrix_t{});
        } else {
          //auto r_mat = ad.cdag_matrix(lin_index, i);
          auto r_mat = (op_dag ? ad.cdag_matrix(lin_index, i) : ad.c_matrix(lin_index, i));
          res.emplace_back(l_bl, l_mat * r_mat);
        }
      }
    }
    return res;

  } // namespace inchworm

  // Calculate the Frobenius norm of the g_frame block diagonal matrix:
  double frobenius_norm(g_frame_t const &g_frame) {
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

  // Calculate the Frobenius norm of the u_partial block diagonal matrix:
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

  void print(frame_t const &u_frame, double factor) {
    for (auto block : u_frame) {
      block *= factor;
      std::cout << block;
    }
    std::cout << "\n";
  }

} // namespace inchworm
