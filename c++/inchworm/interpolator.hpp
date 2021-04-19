#pragma once

#include "./types.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include <nda/array_adapter.hpp>

#include <fmt/core.h>

namespace inchworm {

  class interpolator_t {

    public:
    interpolator_t() = default;

    interpolator_t(u_tau_t const &u_tau, long n_tau)
       : n_blocks(u_tau.size()), n_tau(n_tau), datx(n_tau), daty(n_blocks), interp(n_blocks), acc(n_blocks) {
      EXPECTS(n_tau >= 2);

      for (auto n : range(n_tau)) datx[n] = u_tau[0].mesh()[n];

      for (auto bl : range(n_blocks)) {
        daty[bl] = u_tau[bl].data()(range(n_tau), nda::ellipsis());

        interp[bl] = nda::array<gsl_interp *, 2>{u_tau[bl].target_shape()};
        acc[bl]    = nda::array<gsl_interp_accel *, 2>{u_tau[bl].target_shape()};

        for (auto [i, j] : product_range(u_tau[bl].target_shape())) {
          interp[bl](i, j) = (n_tau == 2) ? gsl_interp_alloc(gsl_interp_linear, n_tau) : gsl_interp_alloc(gsl_interp_cspline, n_tau);
          acc[bl](i, j)    = gsl_interp_accel_alloc();
          gsl_interp_init(interp[bl](i, j), datx.data(), daty[bl](range(), i, j).data(), n_tau);
        }
      }
    }

    // This object holds raw pointers and can only be move-constructed
    interpolator_t(interpolator_t const &) = delete;
    interpolator_t(interpolator_t &&)      = default;

    // This object holds raw pointers and can only be move-assigned
    interpolator_t &operator=(interpolator_t const &) = delete;
    interpolator_t &operator=(interpolator_t &&) = default;

    ~interpolator_t() {
      for (auto bl : range(n_blocks)) {
        if (not interp.empty())
          for (auto *ptr : interp[bl]) gsl_interp_free(ptr);
        if (not acc.empty())
          for (auto *ptr : acc[bl]) gsl_interp_accel_free(ptr);
      }
    }

    scalar_t operator()(int bl, double tau, int i, int j) const {
      EXPECTS(0 <= tau && tau <= datx[n_tau - 1]);
      double res = gsl_interp_eval(interp[bl](i, j), datx.data(), daty[bl](range(), i, j).data(), tau, acc[bl](i, j));
      if (interpolation_failed) { // Store data to file and abort
        {
          auto f = h5::file("interp_debug.h5", 'w');
          h5::write(f, "xvals", datx);
          h5::write(f, "yvals", daty[bl](range(), i, j));
          h5::write(f, "tau", tau);
        }
        std::abort();
      };
      return res;
    }

    matrix_t operator()(int bl, double tau) const {
      EXPECTS(0 <= tau && tau <= datx[n_tau - 1]);
      return nda::array_adapter{interp[bl].shape(), [&](int i, int j) { return (*this)(bl, tau, i, j); }};
    }

    frame_t operator()(double tau) const {
      EXPECTS(0 <= tau && tau <= datx[n_tau - 1]);
      return nda::array_adapter{std::array{n_blocks}, [&](int bl) { return (*this)(bl, tau); }};
    }

    private:
    int n_blocks = 0;
    int n_tau    = 0;

    nda::array<double, 1> datx;
    nda::array<nda::array<scalar_t, 3, nda::F_layout>, 1> daty;

    nda::array<nda::array<gsl_interp *, 2>, 1> interp;
    nda::array<nda::array<gsl_interp_accel *, 2>, 1> acc;

    inline static bool interpolation_failed = false;
    inline static auto custom_error_handler = [](const char *reason, const char *file, int line, int gsl_errno) {
      fmt::print("Interpolation failed (ErrNo: {}, File: {}:{})\n", gsl_errno, file, line);
      interpolation_failed = true;
    };
    inline static auto default_error_handler = gsl_set_error_handler(custom_error_handler);
  };

} // namespace inchworm
