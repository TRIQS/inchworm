#include "./u_frame.hpp"
#include "./util.hpp"

namespace inchworm {

  // --------------- General frame / u_partial functionality ---------------

  frame_t make_frame(std::vector<int> const &shape_of_frame) {
    auto res = frame_t{shape_of_frame.size()};

    for (auto [bl, n] : enumerate(shape_of_frame)) {
      res[bl] = matrix_t{n, n};
      res[bl] = 0.;
    }
    return res;
  }

  frame_t make_frame(gf_struct_t const &gf_struct) {
    auto res = frame_t{gf_struct.size()};

    for (auto bl : range(gf_struct.size())) {
      auto &[blname, blsize] = gf_struct[bl];
      res[bl]                = matrix_t{blsize, blsize};
      res[bl]                = 0.;
    }
    return res;
  }

  frame_t make_frame(u_partial_t const &up) {
    frame_t res{up.size()};
    for (auto bl : range(up.size())) {
      auto &[tbl, mat] = up[bl];
      EXPECTS(bl == tbl || tbl == -1);
      res[bl] = mat;
    }
    return res;
  }

  u_partial_t make_u_partial(frame_t const &u) {
    u_partial_t res(u.size());
    for (int bl = 0; bl < u.size(); ++bl) { res[bl] = {bl, u[bl]}; }
    return res;
  }

  u_partial_t operator*(u_partial_t const &l, u_partial_t const &r) {

    EXPECTS(l.size() == r.size());

    auto res = u_partial_t(l.size());
    for (int i = 0; i < l.size(); ++i) {
      // u[bl] = up[bl, blp] * up[blp, bl]
      // l_bl <- r_bl <- i
      auto &[r_bl, r_mat] = r[i];
      if (r_bl == -1)
        res[i] = {-1, matrix_t{}};
      else {
        auto &[l_bl, l_mat] = l[r_bl];
        if (l_bl == -1)
          res[i] = {-1, matrix_t{}};
        else
          res[i] = {l_bl, l_mat * r_mat};
      }
    }
    return res;
  }

  double frobenius_norm(frame_t const &frame) {
    double val = 0;
    for (auto const &mat : frame) {
      double norm = frobenius_norm(mat);
      val += norm * norm;
    }
    return std::sqrt(val);
  }

  // calculate the trace of the u_frame block diagonal matrix:
  double trace(frame_t const &frame) {
    double val = 0;
    for (auto const &B : frame) val += trace(B);
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

  frame_t make_g_frame_from_l_and_r(atom_diag const &ad_imp, gf_struct_t const &gf_struct, u_partial_t const &l, u_partial_t const &r) {

    frame_t g_frame = make_frame(gf_struct);

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

  // --------------- Atom Diag specific functions -----------------------

  frame_t make_bare_u_frame(atom_diag const &ad, double tau, bool set_gs_to_0) {
    auto u_frame = make_frame(ad.get_subspace_dims());
    u_frame[bl_] << 0.;
    for (auto [bl, bl_size] : enumerate(ad.get_subspace_dims()))
      for (int i : range(bl_size)) u_frame[bl](i, i) = std::exp(-tau * (ad.get_eigenvalue(bl, i) + (set_gs_to_0 ? 0. : ad.get_gs_energy())));
    return u_frame;
  }

  u_tau_t make_ED_propagator(atom_diag const &ad_tot, atom_diag const &ad_atom, atom_diag const &ad_bath, double beta, int n_tau) {

    auto u_tau = u_tau_t{{beta, Fermion, n_tau}, ad_atom.get_subspace_dims()};

    for (int i_tau = 0; i_tau < n_tau; i_tau++) {
      double dtau  = beta * i_tau / (n_tau - 1.);
      auto Z_bath  = trace(ad_bath, [beta](double E) { return std::exp(-beta * E); });
      auto u_frame = partial_trace_bath(ad_tot, ad_atom, ad_bath, beta, dtau) / Z_bath;
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

  // --------------- Block Gf specific functions -----------------------

  frame_t get_frame(u_tau_t const &u_tau, int idx) {
    frame_t res{u_tau.size()};
    for (auto bl : range(u_tau.size())) res[bl] = u_tau[bl][idx];
    return res;
  }

  double relative_distance(u_tau_t const &l, u_tau_t const &r) {
    double dist = 0.0;
    for (int i = 0; i < l[0].mesh().size(); ++i) dist = std::max(dist, relative_distance(get_frame(l, i), get_frame(r, i)));
    return dist;
  }

} // namespace inchworm
