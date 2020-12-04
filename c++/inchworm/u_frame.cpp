#include "./u_frame.hpp"
#include "./util.hpp"

namespace inchworm {

  // --------------- General frame / u_partial functionality ---------------

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

  frame_t make_zero_frame(std::vector<int> const &shape_of_frame) {
    auto res = frame_t{shape_of_frame.size()};

    for (auto [bl, n] : enumerate(shape_of_frame)) {
      res[bl] = matrix_t::zeros({n, n});
    }
    return res;
  }

  frame_t make_zero_frame(gf_struct_t const &gf_struct) {
    auto res = frame_t{gf_struct.size()};

    for (auto bl : range(gf_struct.size())) {
      auto &[blname, blsize] = gf_struct[bl];
      res[bl]                = matrix_t::zeros({blsize, blsize});
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
