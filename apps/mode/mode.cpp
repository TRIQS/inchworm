#include "./mode.hpp"

void base_mode::construct_Hubbard() {
  // prepare input
  std::tie(vi, wi)                          = select_quadrature_GK(n_GK, 0, 1);
  std::tie(Delta_tau, ad_imp, u_tau, G_tau) = test_setup(n_site, n_bath, n_spin, U, mu, t, cp, theta, epsilon);
  u_interpolator                            = interpolator_t<scalar_t>(u_tau, u_tau[0].mesh().size());
  u_tau_max_zeroth_order                    = u_interpolator(tau_max - tau_split) * u_interpolator(tau_split); //oder 0 result
  n_bl                                      = cp.gf_struct.size();
  all_d_ops.resize(n_bl, {});
  all_d_dag_ops.resize(n_bl, {});
  for (auto [bl, bl_pair] : enumerate(cp.gf_struct)) {
    auto [bl_name, bl_size] = bl_pair;
    block_shape.push_back(bl_size);
  }
  fops = fundamental_operator_set{cp.gf_struct};
  for (auto [bl, bl_pair] : enumerate(cp.gf_struct)) {
    auto [bl_name, bl_size] = bl_pair;
    all_d_ops[bl].clear();
    all_d_dag_ops[bl].clear();
    for (auto idx : range(bl_size)) {
      all_d_ops[bl].emplace_back(0.0, false, fops[{bl_name, idx}], bl, idx);
      all_d_dag_ops[bl].emplace_back(0.0, true, fops[{bl_name, idx}], bl, idx);
    }
  }
  if (debug) {
    std::cout << "Delta_tau shape:" << std::endl;
    print_block_shape(Delta_tau);
    std::cout << "G_tau shape:" << std::endl;
    print_block_shape(G_tau);
    std::cout << "u_tau shape:" << std::endl;
    print_block_shape(u_tau);
  }
}

void base_mode::print_summary() {
  int i = subspace_index / u_tau[bl_index].target_shape()[0];
  int j = subspace_index % u_tau[bl_index].target_shape()[0];

  std::cout << "u_tau_max exact: " << std::setw(10) << u_interpolator(tau_max)[bl_index](i, j) << std::endl;
  std::cout << std::left << std::setw(10) << "order" << std::setw(30) << "value" << std::setw(30) << "time(s)" << std::endl;
  std::cout << std::setw(10) << "0" << std::setw(30) << u_tau_max_zeroth_order[bl_index](i, j) << std::endl;
  for (int i = 0; i < order_list.size(); i++) {
    std::cout << std::setw(10) << order_list[i] << std::setw(30) << integral_order_list[i] << std::setw(30) << calculation_time_list[i] << std::endl;
  }
  double sum_value = u_tau_max_zeroth_order[bl_index](i, j) + std::accumulate(integral_order_list.begin(), integral_order_list.end(), 0.0);
  double sum_time  = std::accumulate(calculation_time_list.begin(), calculation_time_list.end(), 0.0);
  std::cout << std::setw(10) << "sum:" << std::setw(30) << sum_value << std::setw(10) << sum_time << std::endl;
}