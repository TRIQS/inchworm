#pragma once
#include <h5/h5.hpp>

inline void h5_save_params(const ModeBase *mode, h5::group h5group, std::string subgroup_name) {
  auto grp = h5group.create_group(subgroup_name);
  h5_write(grp, "n_tot", mode->sp.n_tot);
  h5_write(grp, "n_tau_linear", mode->sp.n_tau_linear);
  h5_write(grp, "order_Chebyshev", mode->sp.order_Chebyshev);
}

inline void h5_save_propagator(const ModeBase *mode, h5::group h5group, std::string subgroup_name) {
  auto grp = h5group.create_group(subgroup_name);
  h5_write(grp, "tau_grid", mode->sp.grid);
  for (auto bl : range(mode->mp.ad_imp.n_subspaces())) {
    for (auto i : range(mode->mp.ad_imp.get_subspace_dim(bl))) {
      for (auto j : range(mode->mp.ad_imp.get_subspace_dim(bl))) {
        std::string dataset_name = fmt::format("u_tau_{}_{}_{}", bl, i, j);
        std::vector<double> data;
        for (size_t i_tau = 0; i_tau < mode->sp.grid.size(); i_tau++) { data.push_back(mode->sr.u_tau[bl][i_tau](i, j)); }
        h5_write(grp, dataset_name, data);
      }
    }
  }
}

inline void h5_save_cheb_coeff(const ModeBase *mode, h5::group h5group, std::string subgroup_name) {
  auto grp                    = h5group.create_group(subgroup_name);
  auto chebyshov_coefficients = mode->sr.u_interpolator.get_cheb_coeffs();
  for (auto bl : range(mode->mp.ad_imp.n_subspaces())) {
    for (auto i : range(mode->mp.ad_imp.get_subspace_dim(bl))) {
      for (auto j : range(mode->mp.ad_imp.get_subspace_dim(bl))) {
        for (auto tau_interval : range(mode->sp.n_tau_linear - 1)) {
          std::cout << "chebyshov_coefficients[" << bl << "][" << tau_interval << "](" << i << "," << j << ")" << std::endl;
          auto coeffs              = chebyshov_coefficients[bl][tau_interval](i, j);
          std::string dataset_name = fmt::format("cheb_coeff_{}_{}_{}_{}", bl, i, j, tau_interval);
          h5_write(grp, dataset_name, coeffs);
        }
      }
    }
  }
}
