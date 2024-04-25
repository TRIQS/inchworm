#pragma once

#include "./types.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_chebyshev.h>

#include <nda/array_adapter.hpp>

#include <fmt/core.h>

#include "./atom_diag.hpp"

using namespace std::complex_literals;

namespace inchworm {

  template <typename T> class interpolator_cspline_t;

  template <> class interpolator_cspline_t<double> {

    public:
    interpolator_cspline_t() = default;

    interpolator_cspline_t(u_tau_t::real_t const &u_tau, long n_tau)
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

    interpolator_cspline_t(u_tau_t::real_t const &u_tau, std::vector<double> const &_datx)
       : n_blocks(u_tau.size()), datx(_datx.size()), daty(n_blocks), interp(n_blocks), accel_ptr(gsl_interp_accel_alloc()) {
      long n_tau = datx.size();
      EXPECTS(n_tau >= 2);
      for(auto n: range(n_tau)) datx[n] = _datx[n];
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
    interpolator_cspline_t(interpolator_cspline_t const &) = delete;
    interpolator_cspline_t(interpolator_cspline_t &&)      = default;

    // This object holds raw pointers and can only be move-assigned
    interpolator_cspline_t &operator=(interpolator_cspline_t const &) = delete;
    interpolator_cspline_t &operator=(interpolator_cspline_t &&)      = default;

    ~interpolator_cspline_t() {
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

  template <> class interpolator_cspline_t<dcomplex> {

    public:
    interpolator_cspline_t() = default;

    interpolator_cspline_t(u_tau_t const &u_tau, long n_tau)
       : n_blocks(u_tau.size()), interpolator_real(real(u_tau), n_tau), interpolator_imag(imag(u_tau), n_tau) {}

    // This object can only be move-constructed as it has members that hold raw pointers
    interpolator_cspline_t(interpolator_cspline_t const &) = delete;
    interpolator_cspline_t(interpolator_cspline_t &&)      = default;

    // This object can only be move-assigned as it has members that hold raw pointers
    interpolator_cspline_t &operator=(interpolator_cspline_t const &) = delete;
    interpolator_cspline_t &operator=(interpolator_cspline_t &&)      = default;

    dcomplex operator()(int bl, double tau, int i, int j) const { return {interpolator_real(bl, tau, i, j), interpolator_imag(bl, tau, i, j)}; }

    nda::matrix<dcomplex> operator()(int bl, double tau) const { return interpolator_real(bl, tau) + 1i * interpolator_imag(bl, tau); }

    nda::array<nda::matrix<dcomplex>, 1> operator()(double tau) const {
      return nda::array_adapter{std::array{n_blocks}, [&](int bl) { return (*this)(bl, tau); }};
    }

    private:
    int n_blocks = 0;
    interpolator_cspline_t<double> interpolator_real;
    interpolator_cspline_t<double> interpolator_imag;
  };

  struct my_f_params {
    long bl;
    long i;
    long j;
    interpolator_cspline_t<double> *interp_cspline;
  };
  inline double my_f(double x, void *p) {
    my_f_params *params = (my_f_params *)p;
    std::cout << "x = " << x << std::endl;

    return params->interp_cspline->operator()(params->bl, x, params->i, params->j);
  }
  template <typename T> class interpolator_linear_Chebyshev_t;

  template <> class interpolator_linear_Chebyshev_t<double> {
    public:
    interpolator_linear_Chebyshev_t() = default;
    // n_tau_linear refers to the number of linear grid points, therefore n_tau_linear-1 linear intervals.
    // order_Chebyshev refers to the Chebyshev order within each of the linear interval, i.e., order_Chebyshev+1 points is needed within the interval.
    // n_tot refers to the total number of grid points. Therefore, each linear interval has (n_tot-n_tau_linear)/(n_tau_linear-1) cbspline points.
    interpolator_linear_Chebyshev_t(u_tau_t::real_t const &u_tau, long n_tot, long n_tau_linear, int order_Chebyshev)
       : n_blocks(u_tau.size()),
         n_tot(n_tot),
         n_tau_linear(n_tau_linear),
         order_Chebyshev(order_Chebyshev),
         datx_linear(n_tau_linear),
         daty(n_blocks),
         //  datx_Chebyshev(n_tau - 1, order + 1),
         interp(n_blocks) {
      EXPECTS(n_tau_linear >= 2);
      EXPECTS(order_Chebyshev >= 2);
      EXPECTS(n_tot >= n_tau_linear);
      // UGLY WORK AROUND:
      // We assume when order_Chebyshev is non-zero, i.e., when we use the linear-Chebyshev grid, the grid points in u_tau are on the linear-Chebyshev grid. It will not match the points in the u_tau[0].mesh() which contains a uniform grid.
      EXPECTS(n_tot == static_cast<long>(n_tau_linear+(n_tau_linear-1)*(order_Chebyshev+1)));
      n_inter = static_cast<long>((n_tot - n_tau_linear) / (n_tau_linear - 1));
      for (auto bl : range(n_blocks)) { daty[bl] = u_tau[bl].data()(range(n_tot), nda::ellipsis()); }
      for (auto n : range(n_tau_linear)) {
        datx_linear[n] = u_tau[0].mesh()[n * (n_inter + 1)];
        std::cout << "datx_linear[" << n << "] = " << datx_linear[n] << std::endl;
      }
      // if (datx_linear[n_tau_linear - 1] != u_tau[0].mesh()[n_tot - 1]) {
      //   std::cout << "datx_linear[" << n_tau_linear - 1 << "] = " << datx_linear[n_tau_linear - 1] << " != " << u_tau[0].mesh()[n_tot - 1]
      //             << std::endl;
      //   std::cout << "updated the last grid point" << std::endl;
      //   datx_linear[n_tau_linear - 1] = u_tau[0].mesh()[n_tot - 1];
      // }
      for (auto bl : range(n_blocks)) {
        interp[bl] = nda::array<nda::array<gsl_cheb_series *, 2>, 1>(n_tau_linear - 1);
        for (auto n : range(n_tau_linear - 1)) { interp[bl][n] = nda::array<gsl_cheb_series *, 2>{u_tau[bl].target_shape()}; }
      }
      // UGLY WORK AROUND
      std::vector<double> grid = generate_linear_Chebyshev_grid(0,datx_linear[n_tau_linear-1], n_tau_linear, order_Chebyshev);
      std::cout << "linear-Chebyshev grid:"<< std::endl;
      for(auto x: grid) std::cout << x << std::endl;
      interpolator_cspline_t<double> interp_cspline(u_tau, grid);
      for (auto n : range(n_tau_linear - 1)) {
        double a = datx_linear[n];
        double b = datx_linear[n + 1];
        if (a > b) throw std::runtime_error("The grid is not consistent with the linear_Chebyshev interpolation");
        std::cout << "a = " << a << ", b = " << b << std::endl;
        // for (auto i : range(order + 1)) { datx_Chebyshev(n, i) = (a + b) / 2 + (b - a) / 2 * std::cos(M_PI * (2 * i + 1) / (2 * (order + 1))); }
        for (auto bl : range(n_blocks)) {
          for (auto [i, j] : product_range(u_tau[bl].target_shape())) {
            // std::cout << "bl = " << bl << ", n = " << n << ", i = " << i << ", j = " << j << std::endl;
            my_f_params params;
            params.bl             = bl;
            params.i              = i;
            params.j              = j;
            params.interp_cspline = &interp_cspline;
            gsl_function F;
            F.function          = &my_f;
            F.params            = &params;
            interp[bl][n](i, j) = gsl_cheb_alloc(order_Chebyshev);
            gsl_cheb_init(interp[bl][n](i, j), &F, a, b);
            // std::cout << "cspline a: "<< params.interp_cspline(bl, a, i, j) << std::endl;
            // std::cout << "cspline b: "<< params.interp_cspline(bl, b, i, j) << std::endl;
            // std::cout << "chebyshev a: "<< (*this)(bl, a, i, j) << std::endl;
            // std::cout << "chebyshev b: "<< (*this)(bl, b, i, j) << std::endl;
          }
        }
      }
      std::cout << "interpolator_linear_Chebyshev_t is constructed" << std::endl;
    }

    // This object holds raw pointers and can only be move-constructed
    interpolator_linear_Chebyshev_t(interpolator_linear_Chebyshev_t const &) = delete;
    interpolator_linear_Chebyshev_t(interpolator_linear_Chebyshev_t &&)      = default;

    // This object holds raw pointers and can only be move-assigned
    interpolator_linear_Chebyshev_t &operator=(interpolator_linear_Chebyshev_t const &) = delete;
    interpolator_linear_Chebyshev_t &operator=(interpolator_linear_Chebyshev_t &&)      = default;

    ~interpolator_linear_Chebyshev_t() {
      for (auto bl : range(n_blocks)) {
        for (auto n : range(n_tau_linear - 1)) {
          if (not interp.empty()) {
            for (auto *ptr : interp[bl][n]) {
              gsl_cheb_free(ptr);
            }
          }
        }
      }
    }

    double operator()(int bl, double tau, int i, int j) const {
      EXPECTS(0 <= tau && tau <= datx_linear[datx_linear.size() - 1]);
      // find the interval
      int n = 0;
      while (n < n_tau_linear - 1 && datx_linear[n + 1] < tau) { n++; }
      if (datx_linear[n] == tau) { return daty[bl](n * (n_inter + 1), i, j); }
      if (datx_linear[datx_linear.size() - 1] == tau) { return daty[bl](n_tot - 1, i, j); }
      if (datx_linear[n + 1] == tau) { return daty[bl]((n + 1) * (n_inter + 1), i, j); }
      return gsl_cheb_eval(interp[bl][n](i, j), tau);
    }

    nda::matrix<double> operator()(int bl, double tau) const {
      EXPECTS(0 <= tau && tau <= datx_linear[datx_linear.size() - 1]);
      return nda::array_adapter{interp[bl][0].shape(), [&](int i, int j) { return (*this)(bl, tau, i, j); }};
    }

    nda::array<nda::matrix<double>, 1> operator()(double tau) const {
      EXPECTS(0 <= tau && tau <= datx_linear[datx_linear.size() - 1]);
      return nda::array_adapter{std::array{n_blocks}, [&](int bl) { return (*this)(bl, tau); }};
    }

    private:
    int n_blocks        = 0;
    long n_tot          = 0; // number of total grid points
    long n_tau_linear   = 0; // number of linear grid points
    long n_inter        = 0; // number of cspine points within each linear interval
    int order_Chebyshev = 0; // Chebyshev order
    nda::array<double, 1> datx_linear;
    nda::array<nda::array<double, 3, nda::F_layout>, 1> daty;
    // nda::array
    // nda::array<double, 2> datx_Chebyshev;
    // block_index, linear_grid_index, block_sub_index
    nda::array<nda::array<nda::array<gsl_cheb_series *, 2>, 1>, 1> interp;
  };

  template <> class interpolator_linear_Chebyshev_t<dcomplex> {
    public:
    interpolator_linear_Chebyshev_t() = default;
    interpolator_linear_Chebyshev_t(u_tau_t const &u_tau, long n_tot, long n_tau_linear, int order_Chebyshev) : n_blocks(u_tau.size()) {
      std::cout << "interpolator_linear_Chebyshev_t<dcomplex> is started" << std::endl;
      interpolator_real = interpolator_linear_Chebyshev_t<double>(real(u_tau), n_tot, n_tau_linear, order_Chebyshev);
      interpolator_imag = interpolator_linear_Chebyshev_t<double>(imag(u_tau), n_tot, n_tau_linear, order_Chebyshev);
      std::cout << "interpolator_linear_Chebyshev_t<dcomplex> is constructed" << std::endl;
    }

    // This object can only be move-constructed as it has members that hold raw pointers
    interpolator_linear_Chebyshev_t(interpolator_linear_Chebyshev_t const &) = delete;
    interpolator_linear_Chebyshev_t(interpolator_linear_Chebyshev_t &&)      = default;

    // This object can only be move-assigned as it has members that hold raw pointers
    interpolator_linear_Chebyshev_t &operator=(interpolator_linear_Chebyshev_t const &) = delete;
    interpolator_linear_Chebyshev_t &operator=(interpolator_linear_Chebyshev_t &&)      = default;

    dcomplex operator()(int bl, double tau, int i, int j) const { return {interpolator_real(bl, tau, i, j), interpolator_imag(bl, tau, i, j)}; }

    nda::matrix<dcomplex> operator()(int bl, double tau) const { return interpolator_real(bl, tau) + 1i * interpolator_imag(bl, tau); }

    nda::array<nda::matrix<dcomplex>, 1> operator()(double tau) const {
      return nda::array_adapter{std::array{n_blocks}, [&](int bl) { return (*this)(bl, tau); }};
    }

    private:
    int n_blocks = 0;
    interpolator_linear_Chebyshev_t<double> interpolator_real;
    interpolator_linear_Chebyshev_t<double> interpolator_imag;
  };
} // namespace inchworm
