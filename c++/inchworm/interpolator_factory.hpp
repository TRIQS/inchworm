#pragma once

#include "./types.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include <nda/array_adapter.hpp>

#include <fmt/core.h>

using namespace std::complex_literals;

namespace inchworm {

  template <typename S> class interpolator_cb_t;

  template <> class interpolator_cb_t<double> {

    public:
    interpolator_cb_t() = default;

    interpolator_cb_t(u_tau_t::real_t const &u_tau, long n_tau)
       : n_blocks(u_tau.size()), datx(n_tau), daty(n_blocks), interp(n_blocks), accel_ptr(gsl_interp_accel_alloc()) {
      EXPECTS(n_tau >= 2);

      for (auto n : range(n_tau)) datx[n] = u_tau[0].mesh()[n];

      for (auto bl : range(n_blocks)) {
        daty[bl]   = u_tau[bl].data()(range(n_tau), nda::ellipsis());
        interp[bl] = nda::array<gsl_interp *, 2>{u_tau[bl].target_shape()};

        for (auto [i, j] : product_range(u_tau[bl].target_shape())) {
          interp[bl](i, j) = (n_tau == 2) ? gsl_interp_alloc(gsl_interp_linear, n_tau) : gsl_interp_alloc(gsl_interp_cspline, n_tau);
          // interp[bl](i, j) = gsl_interp_alloc(gsl_interp_linear, n_tau);
          gsl_interp_init(interp[bl](i, j), datx.data(), daty[bl](range::all, i, j).data(), n_tau);
        }
      }
    }

    // This object holds raw pointers and can only be move-constructed
    interpolator_cb_t(interpolator_cb_t const &) = delete;
    interpolator_cb_t(interpolator_cb_t &&)      = default;

    // This object holds raw pointers and can only be move-assigned
    interpolator_cb_t &operator=(interpolator_cb_t const &) = delete;
    interpolator_cb_t &operator=(interpolator_cb_t &&)      = default;

    ~interpolator_cb_t() {
      for (auto bl : range(n_blocks)) {
        if (not interp.empty())
          for (auto *ptr : interp[bl]) gsl_interp_free(ptr);
      }
    }

    double operator()(int bl, double tau, int i, int j) const {
      EXPECTS(0 <= tau && tau <= datx[datx.size() - 1]);
      double res = gsl_interp_eval(interp[bl](i, j), datx.data(), daty[bl](range::all, i, j).data(), tau, accel_ptr.get());
      if (interpolation_failed) { // Store data to file and abort
        {
          auto f = h5::file("interp_debug.h5", 'w');
          h5::write(f, "xvals", datx);
          h5::write(f, "yvals", daty[bl](range::all, i, j));
          h5::write(f, "tau", tau);
        }
        std::abort();
      };
      return res;
    }

    nda::matrix<double> operator()(int bl, double tau) const {
      EXPECTS(0 <= tau && tau <= datx[datx.size() - 1]);
      return nda::array_adapter{interp[bl].shape(), [&](int i, int j) { return (*this)(bl, tau, i, j); }};
    }

    nda::array<nda::matrix<double>, 1> operator()(double tau) const {
      EXPECTS(0 <= tau && tau <= datx[datx.size() - 1]);
      return nda::array_adapter{std::array{n_blocks}, [&](int bl) { return (*this)(bl, tau); }};
    }

    private:
    int n_blocks = 0;

    nda::array<double, 1> datx;
    nda::array<nda::array<double, 3, nda::F_layout>, 1> daty;

    nda::array<nda::array<gsl_interp *, 2>, 1> interp;

    struct accel_deleter {
      void operator()(auto *p) noexcept { gsl_interp_accel_free(p); };
    };
    std::unique_ptr<gsl_interp_accel, accel_deleter> accel_ptr = {};

    inline static bool interpolation_failed       = false;
    inline static auto const custom_error_handler = [](const char *, const char *file, int line, int gsl_errno) {
      fmt::print("Interpolation failed (ErrNo: {}, File: {}:{})\n", gsl_errno, file, line);
      interpolation_failed = true;
    };
    inline static auto const default_error_handler = gsl_set_error_handler(custom_error_handler);
  };

  template <> class interpolator_cb_t<dcomplex> {

    public:
    interpolator_cb_t() = default;

    interpolator_cb_t(u_tau_t const &u_tau, long n_tau)
       : n_blocks(u_tau.size()), interpolator_real(real(u_tau), n_tau), interpolator_imag(imag(u_tau), n_tau) {}

    // This object can only be move-constructed as it has members that hold raw pointers
    interpolator_cb_t(interpolator_cb_t const &) = delete;
    interpolator_cb_t(interpolator_cb_t &&)      = default;

    // This object can only be move-assigned as it has members that hold raw pointers
    interpolator_cb_t &operator=(interpolator_cb_t const &) = delete;
    interpolator_cb_t &operator=(interpolator_cb_t &&)      = default;

    dcomplex operator()(int bl, double tau, int i, int j) const { return {interpolator_real(bl, tau, i, j), interpolator_imag(bl, tau, i, j)}; }

    nda::matrix<dcomplex> operator()(int bl, double tau) const { return interpolator_real(bl, tau) + 1i * interpolator_imag(bl, tau); }

    nda::array<nda::matrix<dcomplex>, 1> operator()(double tau) const {
      return nda::array_adapter{std::array{n_blocks}, [&](int bl) { return (*this)(bl, tau); }};
    }

    private:
    int n_blocks = 0;
    interpolator_cb_t<double> interpolator_real;
    interpolator_cb_t<double> interpolator_imag;
  };

} // namespace inchworm
