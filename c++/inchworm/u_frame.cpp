#include "./u_frame.hpp"

namespace inchworm {
  //
  u_frame_t make_zero_propagator_frame(atom_diag const &ad) {
    u_frame_t u_frame(ad.n_subspaces());

    for (int bl = 0; bl < ad.n_subspaces(); bl++) {
      u_frame[bl] = matrix_t(ad.get_subspace_dim(bl), ad.get_subspace_dim(bl)); //  use zeros<> ?? check
      u_frame[bl] = 0;
    }
    return u_frame;
  }

  //
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

} // namespace inchworm
