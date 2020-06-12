#include "./u_frame.hpp"

namespace inchworm {
  // Create an empty frame (block diagonal matrix: vector of matrix_t)
  u_frame_t make_zero_propagator_frame(atom_diag const &ad) {
    u_frame_t u_frame(ad.n_subspaces());

    for (int bl = 0; bl < ad.n_subspaces(); bl++) {
      u_frame[bl] = matrix_t(ad.get_subspace_dim(bl), ad.get_subspace_dim(bl)); //  use zeros<> ?? check
      u_frame[bl] = 0;
    }
    return u_frame;
  }

  // initialize bare propagator frame U_0 = exp(-tau H_loc) in the diagonal basis of H_loc
  u_frame_t make_bare_propagator_frame(atom_diag const &ad, double tau, bool set_gs_to_0) {
    u_frame_t u_frame(ad.n_subspaces());

    for (int bl = 0; bl < ad.n_subspaces(); bl++) {
      int dim     = ad.get_subspace_dim(bl);
      u_frame[bl] = matrix_t(dim, dim); //  use zeros<> ?? check
      u_frame[bl] = 0;
      for (int j = 0; j < dim; j++) u_frame[bl](j, j) = std::exp(-tau * (ad.get_eigenvalue(bl, j) + (set_gs_to_0 ? 0. : ad.get_gs_energy())));
    }
    return u_frame;
  }

  u_partial_t make_u_partial(u_frame_t const &u) {
    u_partial_t res;
    for (int i = 0; i < u.size(); ++i) { res.emplace_back(i, u[i]); }
    return res;
  }

  u_frame_t make_u_frame(u_partial_t const &up) {
    u_frame_t res;
    for (int i = 0; i < up.size(); ++i) {
      auto &[bl, mat] = up[i];
      EXPECTS(i == bl);
      EXPECTS(mat.shape()[0] == mat.shape()[1]);
      res.emplace_back(mat);
    }
    return res;
  }

  void multiply_assign(u_frame_t &u_frame, u_partial_t const &l, u_partial_t const &r) {
    EXPECTS(l.size() == r.size());
    EXPECTS(u_frame.size() == r.size());

    for (int i = 0; i < u_frame.size(); ++i) {
      // u[bl] = up[bl, blp] * up[blp, bl]
      // l_bl <- r_bl <- i
      auto &[r_bl, r_mat] = r[i];
      auto &[l_bl, l_mat] = l[r_bl];
      if (r_bl != -1 && l_bl != -1) {
        EXPECTS(l_bl == i);
        EXPECTS(l_mat.shape()[1] == r_mat.shape()[0]);
        EXPECTS(l_mat.shape()[0] == u_frame[i].shape()[0]);
        EXPECTS(r_mat.shape()[1] == u_frame[i].shape()[1]);
        u_frame[i] = l_mat * r_mat;
      }
    }
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

  // calculate the trace of the u_frame block diagonal matrix:
  double trace(u_frame_t const &u_frame) {
    double val = 0;
    for (auto const &B : u_frame) val += trace(B);
    return val;
  }

  void print(u_frame_t const &u_frame, double factor) {
    for (auto block : u_frame) {
      block *= factor;
      std::cout << block;
    }
    std::cout << "\n";
  }

} // namespace inchworm
