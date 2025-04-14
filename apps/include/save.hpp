#pragma once
#include <h5/h5.hpp>

inline void h5_save_params(const ModeBase *mode, h5::group h5group, std::string subgroup_name) {
  auto grp = h5group.create_group(subgroup_name);
  h5_write(grp, "n_tot", mode->sp.n_tot);
  h5_write(grp, "n_tau_linear", mode->sp.n_tau_linear);
  h5_write(grp, "order_Chebyshev", mode->sp.order_Chebyshev);
  h5_write(grp, "propagator_block_shape", mode->mp.ad_imp.get_subspace_dims());
  h5_write(grp, "gf_block_shape", mode->mp.gf_block_shape);
  h5_write(grp, "n_tau_green", mode->cp.n_tau_green);
}

inline void h5_save_propagator(const ModeBase *mode, h5::group h5group, std::string subgroup_name) {
  auto grp = h5group.create_group(subgroup_name);
  h5_write(grp, "tau_grid", mode->sp.grid);
  for (auto bl : range(mode->mp.ad_imp.n_subspaces())) {
    for (auto i : range(mode->mp.ad_imp.get_subspace_dim(bl))) {
      for (auto j : range(mode->mp.ad_imp.get_subspace_dim(bl))) {
        std::string dataset_name = fmt::format("u_tau_{}_{}{}", bl, i, j);
        std::string dataset_name_re = fmt::format("u_tau_re_{}_{}{}", bl, i, j);
        std::vector<double> data;
        std::vector<double> data_re;
        for (size_t i_tau = 0; i_tau < mode->sp.grid.size(); i_tau++) { data.push_back(mode->sr.u_tau[bl][i_tau](i, j) * std::exp(mode->gp.exponent_u * mode->sp.grid[i_tau])); 
        data_re.push_back(mode->sr.u_tau[bl][i_tau](i, j));
        }
        h5_write(grp, dataset_name, data);
        h5_write(grp, dataset_name_re, data_re);
      }
    }
  }
  h5_write(grp, "Z_energy_shift_correction", mode->gp.Z_energy_shift_correction);
  h5_write(grp, "Z_imp_correction", mode->sr.Z_imp_correction);
  h5_write(grp, "exponent_u", mode->gp.exponent_u);
}

inline void h5_save_propagator_ref(const ModeBase *mode, h5::group h5group, std::string subgroup_name) {
  auto grp = h5group.create_group(subgroup_name);
  h5_write(grp, "tau_grid", mode->sp.grid);
  for (auto bl : range(mode->mp.ad_imp.n_subspaces())) {
    for (auto i : range(mode->mp.ad_imp.get_subspace_dim(bl))) {
      for (auto j : range(mode->mp.ad_imp.get_subspace_dim(bl))) {
        std::string dataset_name = fmt::format("u_tau_{}_{}{}", bl, i, j);
        std::vector<double> data;
        for (size_t i_tau = 0; i_tau < mode->sp.grid.size(); i_tau++) { data.push_back(mode->sr.u_tau_ref[bl][i_tau](i, j)); }
        h5_write(grp, dataset_name, data);
      }
    }
  }
  h5_write(grp, "Z_energy_shift_correction", mode->gp.Z_energy_shift_correction);
  h5_write(grp, "Z_imp_correction", mode->sr.Z_imp_correction);
}

inline void h5_save_cheb_coeff(const ModeBase *mode, h5::group h5group, std::string subgroup_name) {
  auto grp                    = h5group.create_group(subgroup_name);
  auto chebyshev_coefficients = mode->sr.u_interpolator.get_cheb_coeffs();
  for (auto bl : range(mode->mp.ad_imp.n_subspaces())) {
    for (auto i : range(mode->mp.ad_imp.get_subspace_dim(bl))) {
      for (auto j : range(mode->mp.ad_imp.get_subspace_dim(bl))) {
        for (auto tau_interval : range(mode->sp.n_tau_linear - 1)) {
          auto coeffs              = chebyshev_coefficients[bl][tau_interval](i, j);
          std::string dataset_name = fmt::format("cheb_coeff_{}_{}_{}{}", bl, i, j, tau_interval);
          h5_write(grp, dataset_name, coeffs);
        }
      }
    }
  }
}

inline void h5_save_cheb_coeff_ref(const ModeBase *mode, h5::group h5group, std::string subgroup_name) {
  auto grp                    = h5group.create_group(subgroup_name);
  auto chebyshev_coefficients = mode->sr.u_interpolator_ref.get_cheb_coeffs();
  for (auto bl : range(mode->mp.ad_imp.n_subspaces())) {
    for (auto i : range(mode->mp.ad_imp.get_subspace_dim(bl))) {
      for (auto j : range(mode->mp.ad_imp.get_subspace_dim(bl))) {
        for (auto tau_interval : range(mode->sp.n_tau_linear - 1)) {
          auto coeffs              = chebyshev_coefficients[bl][tau_interval](i, j);
          std::string dataset_name = fmt::format("cheb_coeff_{}_{}_{}{}", bl, i, j, tau_interval);
          h5_write(grp, dataset_name, coeffs);
        }
      }
    }
  }
}

inline void h5_save_gf(const ModeBase *mode, h5::group h5group, std::string subgroup_name, g_tau_t const &G_tau, std::vector<double> const &tau_grid) {
  auto grp = h5group.create_group(subgroup_name);
  h5_write(grp, "tau_grid", tau_grid);
  for (int bl = 0; bl < mode->mp.gf_block_shape.size(); bl++) {
    for (auto i : range(mode->mp.gf_block_shape[bl])) {
      for (auto j : range(mode->mp.gf_block_shape[bl])) {
        std::string dataset_name = fmt::format("G_tau_{}_{}{}", bl, i, j);
        std::vector<double> data;
        for (size_t i_tau = 0; i_tau < G_tau[0].mesh().size(); i_tau++) { data.push_back(G_tau[bl][i_tau](i, j)); }
        h5_write(grp, dataset_name, data);
      }
    }
  }
}

inline void h5_save_statistics(const ModeBase *mode, h5::group h5group, std::string subgroup_name) {
  auto grp = h5group.create_group(subgroup_name);
  // sr.statistics is a vector of vector, the first one is for different inchworm/bare step and the second one is the for different order
  // generate subgroup for different components of statistics
  auto grp_func_evals = grp.create_group("func_evals");
  auto grp_warning_same_time = grp.create_group("warning_same_time");
  auto grp_warning_tau_split = grp.create_group("warning_tau_split");
  auto grp_warning_tau_max = grp.create_group("warning_tau_max");
  auto grp_max_diff = grp.create_group("max_diff");
  auto grp_max_auxi_height = grp.create_group("max_auxi_height");
  auto grp_max_pivot_error = grp.create_group("max_pivot_error");
  auto grp_u_tau_sum = grp.create_group("u_tau_sum");
  auto grp_time = grp.create_group("time");
  auto grp_nTCI = grp.create_group("nTCI");
  auto grp_integral_max = grp.create_group("integral_max");
  for (size_t i = 0; i < mode->sr.statistics.size(); i++) {
    std::string subgroup_name = fmt::format("inch{}", i);
    h5_write(grp_func_evals, subgroup_name, mode->sr.statistics[i].func_evals_order);
    h5_write(grp_warning_same_time, subgroup_name, mode->sr.statistics[i].warning_same_time_order);
    h5_write(grp_warning_tau_split, subgroup_name, mode->sr.statistics[i].warning_tau_split_order);
    h5_write(grp_warning_tau_max, subgroup_name, mode->sr.statistics[i].warning_tau_max_order);
    h5_write(grp_max_diff, subgroup_name, mode->sr.statistics[i].max_diff_order);
    h5_write(grp_max_auxi_height, subgroup_name, mode->sr.statistics[i].max_auxi_height_order);
    h5_write(grp_max_pivot_error, subgroup_name, mode->sr.statistics[i].max_error_order);
    h5_write(grp_u_tau_sum, subgroup_name, mode->sr.statistics[i].u_tau_sum_order);
    h5_write(grp_time, subgroup_name, mode->sr.statistics[i].time_order);
    h5_write(grp_nTCI, subgroup_name, mode->sr.statistics[i].nTCI_order);
    h5_write(grp_integral_max, subgroup_name, mode->sr.statistics[i].integral_max_order);
}
}