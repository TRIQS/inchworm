#pragma once

#include "./types.hpp"

#include <nda/array_adapter.hpp>

#include <fmt/core.h>

#include <variant>

#include "interpolator_factory.hpp"

using namespace std::complex_literals;

namespace inchworm {
  enum class interpolation_type { cspline, linear_Chebyshev };

  template <typename T> class interpolator_t {

    public:
    interpolator_t() = default;
    // order is 0 for cubic spline, non-zero for linear-Chebyshev
    interpolator_t(u_tau_t const &u_tau, long n_tau, int order, interpolation_type type) : type(type) {
      if (type == interpolation_type::cspline) {
        if (order != 0) throw std::runtime_error("order must be 0 for cubic spline interpolation");
        interpolator = interpolator_cspline_t<T>(u_tau, n_tau);
      } else if (type == interpolation_type::linear_Chebyshev) {
        if (order == 0) throw std::runtime_error("order must be non-zero for linear-Chebyshev interpolation");
        interpolator = interpolator_linear_Chebyshev_t<T>(u_tau, n_tau, order);
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

    private:
    std::variant<interpolator_cspline_t<T>, interpolator_linear_Chebyshev_t<T>> interpolator;
    interpolation_type type;
  };

} // namespace inchworm
