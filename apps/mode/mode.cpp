#include "./mode.hpp"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

void base_mode::prepare_input() {
  std::tie(tp.vi, tp.wi_v)                              = select_quadrature_GK(tp.n_GK, 0, 1);
  std::tie(mp.Delta_tau, mp.ad_imp, sr.u_tau, sr.G_tau) = test_setup(mp.n_site, mp.n_bath, mp.n_spin, mp.U, mp.mu, mp.t, cp, mp.theta, mp.epsilon);
  sr.u_interpolator                                     = interpolator_t<scalar_t>(sr.u_tau, sr.u_tau[0].mesh().size());
  sr.u_tau_max_zeroth_order = sr.u_interpolator(sp.tau_max - sp.tau_split) * sr.u_interpolator(sp.tau_split); //oder 0 result

  int n_bl = cp.gf_struct.size();
  mp.all_d_ops.resize(n_bl, {});
  mp.all_d_dag_ops.resize(n_bl, {});
  for (auto [bl, bl_pair] : enumerate(cp.gf_struct)) {
    auto [bl_name, bl_size] = bl_pair;
    mp.gf_block_shape.push_back(bl_size);
  }
  mp.n_phi = std::accumulate(mp.gf_block_shape.begin(), mp.gf_block_shape.end(), 0);
  mp.fops  = fundamental_operator_set{cp.gf_struct};
  for (auto [bl, bl_pair] : enumerate(cp.gf_struct)) {
    auto [bl_name, bl_size] = bl_pair;
    mp.all_d_ops[bl].clear();
    mp.all_d_dag_ops[bl].clear();
    for (auto idx : range(bl_size)) {
      mp.all_d_ops[bl].emplace_back(0.0, false, mp.fops[{bl_name, idx}], bl, idx);
      mp.all_d_dag_ops[bl].emplace_back(0.0, true, mp.fops[{bl_name, idx}], bl, idx);
    }
  }
  if (sp.debug > 1) {
    std::cout << "Delta_tau shape:" << std::endl;
    print_block_shape(mp.Delta_tau);
    std::cout << "G_tau shape:" << std::endl;
    print_block_shape(sr.G_tau);
    std::cout << "u_tau shape:" << std::endl;
    print_block_shape(sr.u_tau);
  }
}

void base_mode::print_summary() {
  int i = sp.subspace_index / sr.u_tau[sp.bl_index].target_shape()[0];
  int j = sp.subspace_index % sr.u_tau[sp.bl_index].target_shape()[0];

  std::cout << "u_tau_max exact: " << std::setw(10) << sr.u_interpolator(sp.tau_max)[sp.bl_index](i, j) << std::endl;
  std::cout << std::left << std::setw(10) << "order" << std::setw(30) << "value" << std::setw(30) << "time(s)" << std::endl;
  std::cout << std::setw(10) << "0" << std::setw(30) << sr.u_tau_max_zeroth_order[sp.bl_index](i, j) << std::endl;
  for (int i = 0; i < sp.order_list.size(); i++) {
    std::cout << std::setw(10) << sp.order_list[i] << std::setw(30) << sr.integral_order_list[i] << std::setw(30) << sr.calculation_time_list[i]
              << std::endl;
  }
  double sum_value =
     sr.u_tau_max_zeroth_order[sp.bl_index](i, j) + std::accumulate(sr.integral_order_list.begin(), sr.integral_order_list.end(), 0.0);
  double sum_time = std::accumulate(sr.calculation_time_list.begin(), sr.calculation_time_list.end(), 0.0);
  std::cout << std::setw(10) << "sum:" << std::setw(30) << sum_value << std::setw(10) << sum_time << std::endl;
}

void base_mode::read_json_parameters(std::string json_file_path) {

  namespace pt = boost::property_tree;
  pt::ptree root;
  pt::read_json(json_file_path, root);

  // Read construction parameters
  cp.beta        = root.get<double>("cp.beta");
  cp.n_tau_green = root.get<int>("cp.n_tau_green");
  cp.n_tau_inch  = root.get<int>("cp.n_tau_inch");
  cp.n_tau       = root.get<int>("cp.n_tau");
  for (pt::ptree::value_type &g_s : root.get_child("cp.gf_struct")) {
    std::string name = g_s.first;
    int size         = g_s.second.get_value<int>();
    cp.gf_struct.emplace_back(std::make_pair(name, size));
  }

  // Read model parameters
  mp.n_site = root.get<int>("mp.n_site");
  mp.n_bath = root.get<int>("mp.n_bath");
  mp.n_spin = root.get<int>("mp.n_spin");
  mp.U      = root.get<double>("mp.U");
  mp.mu     = root.get<double>("mp.mu");
  mp.t      = root.get<double>("mp.t");

  int size = root.get_child("mp.epsilon").size();
  mp.epsilon.resize(size);
  int i = 0;
  for (pt::ptree::value_type &ep : root.get_child("mp.epsilon")) {
    mp.epsilon[i] = ep.second.get_value<double>();
    i++;
  }

  int sizex = root.get_child("mp.theta").size();
  int sizey = root.get_child("mp.theta").begin()->second.size();
  mp.theta.resize(sizex, sizey);
  i = 0;
  for (pt::ptree::value_type &th : root.get_child("mp.theta")) {
    int j = 0;
    for (pt::ptree::value_type &th_i : th.second) {
      mp.theta(i, j) = th_i.second.get_value<double>();
      j++;
    }
    i++;
  }

  // Read simulation parameters
  sp.tau_max        = root.get<double>("sp.tau_max");
  sp.tau_split      = root.get<double>("sp.tau_split_ratio") * sp.tau_max;
  sp.bl_index       = root.get<int>("sp.bl_index");
  sp.subspace_index = root.get<int>("sp.subspace_index");

  int debug_level = root.get<int>("sp.debug");
  if (debug_level == 0) {
    sp.debug = debug_t::none;
  } else if (debug_level == 1) {
    sp.debug = debug_t::low;
  } else if (debug_level == 2) {
    sp.debug = debug_t::high;
  } else {
    throw std::runtime_error("Invalid debug level");
  }

  size = root.get_child("sp.order_list").size();
  sp.order_list.resize(size);
  i = 0;
  for (pt::ptree::value_type &order : root.get_child("sp.order_list")) {
    sp.order_list[i] = order.second.get_value<int>();
    i++;
  }

  // Read TCI parameters
  tp.n_GK                 = root.get<int>("tp.n_GK");
  tp.tci_prrlu            = root.get<bool>("tp.tci_prrlu");
  tp.bond_dim             = root.get<int>("tp.bond_dim");
  tp.sweep_bound          = root.get<int>("tp.sweep_bound");
  tp.integral_error_bound = root.get<double>("tp.integral_error_bound");
  tp.pivot_error_bound    = root.get<double>("tp.pivot_error_bound");
  tp.auxi_height          = root.get<double>("tp.auxi_height");
  tp.reltol               = root.get<double>("tp.reltol");
  std::cout << "json parameter file read successfully" << std::endl;
}