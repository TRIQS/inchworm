#include "./mode.hpp"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

void ModeBase::read_json_parameters(std::string json_file_path) {
  namespace pt = boost::property_tree;
  pt::ptree root;
  pt::read_json(json_file_path, root);

  // Read global parameters
  gp.target            = root.get<std::string>("gp.target");
  gp.integrand         = root.get<std::string>("gp.integrand");
  gp.integral_variable = root.get<std::string>("gp.integral_variable");
  gp.tci_shape         = root.get<std::string>("gp.tci_shape");
  gp.trick             = root.get<std::string>("gp.trick");
  gp.model_type        = root.get<int>("gp.model_type");

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
  tp.mapping_v            = root.get<int>("tp.mapping_v");
  tp.tci_prrlu            = root.get<bool>("tp.tci_prrlu");
  tp.bond_dim             = root.get<int>("tp.bond_dim");
  tp.sweep_bound          = root.get<int>("tp.sweep_bound");
  tp.auxi_height          = root.get<double>("tp.auxi_height");
  tp.reltol               = root.get<double>("tp.reltol");
  tp.integral_lower_bound = root.get<double>("tp.integral_lower_bound");
  tp.convergence_bound    = root.get<double>("tp.convergence_bound");
  tp.convergence_iter     = root.get<int>("tp.convergence_iter");
  std::cout << "json parameter file read successfully" << std::endl;
} // end of read_json_parameters

void ModeBase::prepare_input() {
  std::tie(tp.v_value, tp.v_weight) = select_quadrature_GK(tp.n_GK, 0, 1);

  if (gp.model_type == 0) { //model_type 0: discrete bath, where exact results (reference) are available

    // input parameters and exact results
    std::tie(mp.Delta_tau, mp.ad_imp, sr.u_tau_ref, sr.G_tau_ref) =
       test_setup(mp.n_site, mp.n_bath, mp.n_spin, mp.U, mp.mu, mp.t, cp, mp.theta, mp.epsilon);
    sr.u_interpolator_ref     = interpolator_t<scalar_t>(sr.u_tau_ref, sr.u_tau_ref[0].mesh().size());
    sr.partition_function_ref = trace(sr.u_interpolator_ref(cp.beta));

    // structure information about Green's function
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
      print_block_shape(sr.G_tau_ref);
      std::cout << "u_tau shape:" << std::endl;
      print_block_shape(sr.u_tau_ref);
    }
  } else if (gp.model_type == 1) { //model_type 1: continuous bath (read from file)
    std::cerr << "not implemented yet" << std::endl;
    std::exit(EXIT_FAILURE);

  } else if (gp.model_type == 2) { //model_type 2: bethe lattice
    std::cerr << "not implemented yet" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  else{
    std::cerr << "invalid model_type" << std::endl;
    std::exit(EXIT_FAILURE);
  }
} // end of prepare_input

void ModeBase::print_summary() {
  // global parameters
  std::cout << "mode_name: " << mode_name << std::endl;
  std::cout << "target: " << gp.target << std::endl;
  std::cout << "integrand: " << gp.integrand << std::endl;
  std::cout << "integral_variable: " << gp.integral_variable << std::endl;
  std::cout << "tci_shape: " << gp.tci_shape << std::endl;
  std::cout << "trick: " << gp.trick << std::endl;
  if (gp.model_type == 0) {
    std::cout << "model_type: discrete bath" << std::endl;
  } else if (gp.model_type == 1) {
    std::cout << "model_type: continuous bath" << std::endl;
  } else if (gp.model_type == 2) {
    std::cout << "model_type: bethe lattice" << std::endl;
  }

  if (mode_name == "debug") {
    // debug mode now only work on evaluating a single element of the propagator
    int i = sp.subspace_index / sr.u_tau[sp.bl_index].target_shape()[0];
    int j = sp.subspace_index % sr.u_tau[sp.bl_index].target_shape()[0];
    std::cout << "element index: "
              << "(" << i << ", " << j << ")" << std::endl;
    std::cout << "u_tau_max exact: " << std::setw(10) << sr.u_interpolator_ref(sp.tau_max)[sp.bl_index](i, j) << std::endl;
    std::cout << std::left << std::setw(10) << "order" << std::setw(30) << "value" << std::setw(30) << "time(s)"
              << "time_find_pivot(s)" << std::setw(30) << "time_pretrain(s)" << std::setw(30) << "time_train(s)" << std::endl;
    std::cout << std::setw(10) << "0" << std::setw(30) << sr.u_tau_zeroth_order_ref[sp.bl_index](i, j) << std::endl;
    for (int i = 0; i < sp.order_list.size(); i++) {
      std::cout << std::setw(10) << sp.order_list[i] << std::setw(30) << sr.integral_list[i] << std::setw(30) << sr.calculation_time_list[i]
                << std::setw(30) << sr.find_pivot_time_list[i] << std::setw(30) << sr.pretrain_time_list[i] << std::setw(30) << sr.train_time_list[i]
                << std::endl;
    }
    double sum_value = sr.u_tau_max_zeroth_order[sp.bl_index](i, j) + std::accumulate(sr.integral_list.begin(), sr.integral_list.end(), 0.0);
    double sum_time  = std::accumulate(sr.calculation_time_list.begin(), sr.calculation_time_list.end(), 0.0);
    double sum_time_find_pivot = std::accumulate(sr.find_pivot_time_list.begin(), sr.find_pivot_time_list.end(), 0.0);
    double sum_time_pretrain   = std::accumulate(sr.pretrain_time_list.begin(), sr.pretrain_time_list.end(), 0.0);
    double sum_time_train      = std::accumulate(sr.train_time_list.begin(), sr.train_time_list.end(), 0.0);
    std::cout << std::setw(10) << "sum:" << std::setw(30) << sum_value << std::setw(10) << sum_time << std::setw(10) << sum_time_find_pivot
              << std::setw(10) << sum_time_pretrain << std::setw(10) << sum_time_train << std::endl;
  } else if (mode_name == "bare") {
    std::cerr << "not implemented yet" << std::endl;
    std::exit(EXIT_FAILURE);
  } else if (mode_name == "inchworm") {
    std::cerr << "not implemented yet" << std::endl;
    std::exit(EXIT_FAILURE);
  } else {
    std::cerr << "invalid mode_name" << std::endl;
    std::exit(EXIT_FAILURE);
  }
} // end of print_summary
