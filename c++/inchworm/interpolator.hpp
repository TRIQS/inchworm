#pragma once

#include "./types.hpp"
#include <nda/array_adapter.hpp>
#include <fmt/core.h>
#include <variant>
#include "interpolator_factory.hpp"
#include "nvtx.hpp"

using namespace std::complex_literals;

namespace inchworm {
  enum class interpolation_type { cspline, linear_Chebyshev };

  template <typename T> class interpolator_t {

    public:
    interpolator_t() = default;
    // order is 0 for cubic spline, non-zero for linear-Chebyshev
    interpolator_t(u_tau_t const &u_tau, long n_tot, long n_tau, int order, interpolation_type type, std::vector<double> const & grid = std::vector<double>()) : type(type) {
      if (type == interpolation_type::cspline) {
        if (order != -1) throw std::runtime_error("order must be -1 for cubic spline interpolation");
        interpolator = interpolator_cspline_t<T>(u_tau, n_tot);
      } else if (type == interpolation_type::linear_Chebyshev) {
        if (order < 0) throw std::runtime_error("order must be positive integer for linear-Chebyshev interpolation");
        if (grid.size() != n_tot) throw std::runtime_error("grid size must be equal to n_tot");
        interpolator = interpolator_linear_Chebyshev_t<T>(u_tau, n_tot, n_tau, order, grid);
      } else {
        throw std::runtime_error("Only cspline and linear-Chebyshev interpolation are supported");
      }
    }

    // This object holds raw pointers and can only be move-constructed
    interpolator_t(interpolator_t const &) = delete;
    interpolator_t(interpolator_t &&)      = default;

    // This object holds raw pointers and can only be move-assigned
    interpolator_t &operator=(interpolator_t const &) = delete;
    interpolator_t &operator=(interpolator_t &&)      = default;

    ~interpolator_t() = default;

    T operator()(int bl, double tau, int i, int j) const {
      NVTX_RANGE("obtain interpolated points", 0);
      std::cout << "interpolator_t::operator() called" << std::endl;
      if (type == interpolation_type::cspline) {
        return std::get<interpolator_cspline_t<T>>(interpolator)(bl, tau, i, j);
      } else if (type == interpolation_type::linear_Chebyshev) {
        return std::get<interpolator_linear_Chebyshev_t<T>>(interpolator)(bl, tau, i, j);
      } else
        throw std::runtime_error("Only cspline and linear-Chebyshev interpolation are supported");
    }

    nda::matrix<T> operator()(int bl, double tau) const {
      if (type == interpolation_type::cspline) {
        return std::get<interpolator_cspline_t<T>>(interpolator)(bl, tau);
      } else if (type == interpolation_type::linear_Chebyshev) {
        return std::get<interpolator_linear_Chebyshev_t<T>>(interpolator)(bl, tau);
      } else
        throw std::runtime_error("Only cspline and linear-Chebyshev interpolation are supported");
    }

    nda::array<nda::matrix<double>, 1> operator()(double tau) const {
      if (type == interpolation_type::cspline) {
        return std::get<interpolator_cspline_t<T>>(interpolator)(tau);
      } else if (type == interpolation_type::linear_Chebyshev) {
        return std::get<interpolator_linear_Chebyshev_t<T>>(interpolator)(tau);
      } else
        throw std::runtime_error("Only cspline and linear-Chebyshev interpolation are supported");
    }

    nda::array<nda::array<nda::array<nda::array<double, 1>, 2>, 1>, 1> get_cheb_coeffs() const {
      if (type == interpolation_type::linear_Chebyshev) {
        return std::get<interpolator_linear_Chebyshev_t<T>>(interpolator).get_cheb_coeffs();
      } else
        throw std::runtime_error("Only linear-Chebyshev interpolation has Chebyshev coefficients");
    }

    private:
    std::variant<interpolator_cspline_t<T>, interpolator_linear_Chebyshev_t<T>> interpolator;
    interpolation_type type;
  };

} // namespace inchworm
