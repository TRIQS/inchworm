#include "./mode.hpp"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <h5/h5.hpp>
#include "../utility.hpp"

void ModeBase::clear_tci_results() {
  sr.integral_list.clear();
  sr.calculation_time_list.clear();
  sr.find_pivot_time_list.clear();
  sr.pretrain_time_list.clear();
  sr.train_time_list.clear();
} // end of clear_results

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
  gp.do_segment        = root.get<bool>("gp.do_segment");
  try {
    gp.output_prefix = root.get<std::string>("gp.output_prefix");
  } catch (const std::exception &e) {
    gp.output_prefix = "sim"; 
  }

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
  mp.n_site        = root.get<int>("mp.n_site");
  mp.n_bath        = root.get<int>("mp.n_bath");
  mp.n_spin        = root.get<int>("mp.n_spin");
  mp.U             = root.get<double>("mp.U");
  mp.mu            = root.get<double>("mp.mu");
  mp.t             = root.get<double>("mp.t");
  mp.n_omega_bethe = root.get<int>("mp.n_omega_bethe");

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
  sp.tau_max         = root.get<double>("sp.tau_max");
  sp.tau_split_ratio = root.get<double>("sp.tau_split_ratio");
  sp.tau_split       = sp.tau_split_ratio * sp.tau_max;
  sp.bl_index        = root.get<int>("sp.bl_index");
  sp.subspace_index  = root.get<int>("sp.subspace_index");

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
  tp.tci_prrlu            = root.get<int>("tp.tci_prrlu");
  tp.bond_dim             = root.get<int>("tp.bond_dim");
  tp.sweep_bound          = root.get<int>("tp.sweep_bound");
  tp.auxi_height          = root.get<double>("tp.auxi_height");
  tp.reltol               = root.get<double>("tp.reltol");
  tp.fullPiv              = root.get<bool>("tp.fullPiv");
  tp.error_type           = root.get<int>("tp.error_type");
  tp.error_eval           = root.get<int>("tp.error_eval");
  tp.convergence_bound    = root.get<double>("tp.convergence_bound");
  tp.convergence_iter     = root.get<int>("tp.convergence_iter");
  tp.integral_lower_bound = root.get<double>("tp.integral_lower_bound");
  std::cout << "json parameter file read successfully" << std::endl;
} // end of read_json_parameters

hyb_tau_t ModeBase::read_hyb_function(std::string hyb_file_path, model_params_t const &mp, constr_params_t const &cp) {
  if (hyb_file_path.empty()) {
    std::cerr << "hyb_file_path is empty" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::vector<int> block_shape{};
  h5::file file{hyb_file_path, 'r'};
  h5::group grp{file};
  h5_read(grp, "bl_structure", block_shape);
  if (block_shape != mp.gf_block_shape) {
    std::cerr << "block shape in hyb_file_path does not match the gf_block_shape in json parameter file" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::vector<double> tau_grid{};
  h5_read(grp, "tau_grid", tau_grid);
  auto Delta_tau = hyb_tau_t{{cp.beta, Fermion, static_cast<long>(tau_grid.size())}, cp.gf_struct};

  auto Delta_tau_grid = Delta_tau[0].mesh();
  for (int i = 0; i < Delta_tau_grid.size(); i++) {
    if (std::abs(Delta_tau_grid[i] - tau_grid[i]) > 1e-10) {
      std::cerr << "tau_grid in hyb_file_path does not match the tau_grid in json parameter file" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  for (int block = 0; block < cp.gf_struct.size(); block++) {
    for (auto [i, j] : product_range(mp.n_site, mp.n_site)) {
      std::vector<double> Delta_tau_ij{};
      Delta_tau_ij.reserve(tau_grid.size());
      h5_read(grp, "data/" + std::to_string(block) + '_' + std::to_string(i) + std::to_string(j), Delta_tau_ij);
      for (int k = 0; k < Delta_tau_grid.size(); k++) { Delta_tau[block][Delta_tau_grid[k]](i, j) = Delta_tau_ij[k]; }
    }
  }

  return Delta_tau;
}

void ModeBase::prepare_input(std::string hyb_file_path) {
  std::tie(tp.v_value, tp.v_weight) = select_quadrature_GK(tp.n_GK, 0, 1);
  // int n_grid_tanh_sinh  = 50;
  // double h_tanh_sinh    = 4.0 / n_grid_tanh_sinh;
  // auto [xi_old, wi_old] = tanh_sinh_quadrature(0, 1, h_tanh_sinh, n_grid_tanh_sinh);
  // //find all 0 and 1 in xi and remove the corresponding xi and wi
  // std::vector<double> xi, wi;
  // for (int i = 0; i < xi_old.size(); i++) {
  //   if (xi_old[i] != 0 && xi_old[i] != 1) {
  //     tp.v_value.push_back(xi_old[i]);
  //     tp.v_weight.push_back(wi_old[i]);
  //   }
  // }
  // // uniform grid
  // int N = 15;
  // std::vector<double> xi(N), wi(N);
  // for (int i = 0; i < N; i++) {
  //   xi[i] = (i+1) / (N + 1.0);
  //   wi[i] = 1.0 / (N + 1.0);
  // }
  // for(int i = 0; i < N; i++) {
  //   tp.v_value.push_back(xi[i]);
  //   tp.v_weight.push_back(wi[i]);
  // }
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

  sp.order_Chebyshev = 20;
  sp.n_tau_linear    = cp.n_tau_inch; // we temporarily set n_tau_linear to be equal to n_tau_inch
  if (sp.order_Chebyshev != 0) {
    sp.interp_type = interpolation_type::linear_Chebyshev;
    sp.n_tot       = sp.n_tau_linear + (sp.n_tau_linear - 1) * (sp.order_Chebyshev + 1);
  } else {
    sp.interp_type = interpolation_type::cspline;
    sp.n_tot       = sp.n_tau_linear;
  }

  if (gp.model_type == 0) { //model_type 0: discrete bath, where exact results (reference) are available
    // sp.grid = generate_inchworm_grid(0, cp.beta, sp.n_tau_linear, sp.order_Chebyshev);
    // // input parameters and exact results
    // sp.grid_linear = generate_inchworm_grid(0, cp.beta, sp.n_tau_linear, 0);
    std::tie(sr.Z_bath_correction, sr.Z_imp_correction, sr.Z_bath, mp.Delta_tau, mp.ad_imp, sr.u_tau_ref, sr.G_tau_ref) =
       discrete_setup(mp.n_site, mp.n_bath, mp.n_spin, mp.U, mp.mu, mp.t, cp, mp.theta, mp.epsilon, sp.n_tot, sp.n_tau_linear, sp.order_Chebyshev,
                      sp.grid_linear, sp.grid);

    sr.u_interpolator_ref                  = interpolator_t<scalar_t>(sr.u_tau_ref, sp.n_tot, sp.n_tau_linear, sp.order_Chebyshev, sp.interp_type);
    sr.partition_function_ref              = trace(sr.u_interpolator_ref(cp.beta));
    sr.u_tau_zeroth_order_ref              = sr.u_interpolator_ref(sp.tau_max - sp.tau_split) * sr.u_interpolator_ref(sp.tau_split); //oder 0 result
    sr.u_tau_zeroth_order_bare             = make_bare_u_frame(mp.ad_imp, cp.beta);
    sr.partition_function_zeroth_order_ref = trace(sr.u_tau_zeroth_order_bare);
    if (sp.debug > 1) {
      std::cout << "Delta_tau shape:" << std::endl;
      print_block_shape(mp.Delta_tau);
      std::cout << "G_tau shape:" << std::endl;
      print_block_shape(sr.G_tau_ref);
      std::cout << "u_tau shape:" << std::endl;
      print_block_shape(sr.u_tau_ref);
      std::cout << "ad_imp shape:" << std::endl;
      for (auto bl : range(mp.ad_imp.n_subspaces())) { std::cout << "bl: " << bl << ", dim: " << mp.ad_imp.get_subspace_dim(bl) << std::endl; }
    }
  } else if (gp.model_type == 1) { //model_type 1: read hybridization function from input file

    mp.ad_imp                              = imp_setup(mp.n_site, mp.n_spin, mp.U, mp.mu, mp.t, cp);
    mp.Delta_tau                           = read_hyb_function(hyb_file_path, mp, cp);
    sr.u_tau_zeroth_order_bare             = make_bare_u_frame(mp.ad_imp, cp.beta);
    sr.partition_function_zeroth_order_ref = trace(sr.u_tau_zeroth_order_bare);
    sr.partition_function_ref              = 0.0;
    std::tie(sp.grid_linear, sp.grid)      = generate_inchworm_grid(0, cp.beta, sp.n_tau_linear, sp.order_Chebyshev);

  } else if (gp.model_type == 2) { //model_type 2: bethe lattice
    std::tie(mp.Delta_tau, mp.ad_imp)      = bethe_setup(mp.n_site, mp.n_spin, mp.U, mp.mu, mp.t, cp, mp.theta, mp.n_omega_bethe);
    sr.u_tau_zeroth_order_bare             = make_bare_u_frame(mp.ad_imp, cp.beta);
    sr.partition_function_zeroth_order_ref = trace(sr.u_tau_zeroth_order_bare);
    sr.partition_function_ref              = 0.0;
  } else {
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
    int i = sp.subspace_index / sr.u_tau_ref[sp.bl_index].target_shape()[0];
    int j = sp.subspace_index % sr.u_tau_ref[sp.bl_index].target_shape()[0];
    std::cout << "element index: " << "(" << i << ", " << j << ")" << std::endl;
    std::cout << "u_tau_max exact: " << std::setw(10) << sr.u_interpolator_ref(sp.tau_max)[sp.bl_index](i, j) << std::endl;
    std::cout << "u_tau_max exact * Z_imp_correction: " << std::setw(10) << sr.u_interpolator_ref(sp.tau_max)[sp.bl_index](i, j) * sr.Z_imp_correction
              << std::endl;
    std::cout << "u_tau_max exact * Z_bath: " << std::setw(10) << sr.u_interpolator_ref(sp.tau_max)[sp.bl_index](i, j) * sr.Z_bath << std::endl;
    std::cout << "u_tau_max exact * Z_bath * Z_imp_correction*Z_bath_correction: " << std::setw(10)
              << sr.u_interpolator_ref(sp.tau_max)[sp.bl_index](i, j) * sr.Z_bath * sr.Z_imp_correction * sr.Z_bath_correction << std::endl;
    std::cout << std::left << std::setw(10) << "order" << std::setw(30) << "value" << std::setw(30) << "time(s)" << std::setw(30)
              << "time_find_pivot(s)" << std::setw(30) << "time_pretrain(s)" << std::setw(30) << "time_train(s)" << std::endl;
    std::cout << std::setw(10) << "0" << std::setw(30) << sr.u_tau_zeroth_order_ref[sp.bl_index](i, j) << std::endl;
    for (int i = 0; i < sp.order_list.size(); i++) {
      std::cout << std::setw(10) << sp.order_list[i] << std::setw(30) << sr.integral_list[i] << std::setw(30) << sr.calculation_time_list[i]
                << std::setw(30) << sr.find_pivot_time_list[i] << std::setw(30) << sr.pretrain_time_list[i] << std::setw(30) << sr.train_time_list[i]
                << std::endl;
    }
    double sum_value = sr.u_tau_zeroth_order_ref[sp.bl_index](i, j) + std::accumulate(sr.integral_list.begin(), sr.integral_list.end(), 0.0);
    double sum_time  = std::accumulate(sr.calculation_time_list.begin(), sr.calculation_time_list.end(), 0.0);
    double sum_time_find_pivot = std::accumulate(sr.find_pivot_time_list.begin(), sr.find_pivot_time_list.end(), 0.0);
    double sum_time_pretrain   = std::accumulate(sr.pretrain_time_list.begin(), sr.pretrain_time_list.end(), 0.0);
    double sum_time_train      = std::accumulate(sr.train_time_list.begin(), sr.train_time_list.end(), 0.0);
    std::cout << std::setw(10) << "sum:" << std::setw(30) << sum_value << std::setw(30) << sum_time << std::setw(30) << sum_time_find_pivot
              << std::setw(30) << sum_time_pretrain << std::setw(30) << sum_time_train << std::endl;
    std::cout << "sum * Z_imp_correction: " << std::setw(10) << sum_value * sr.Z_imp_correction << std::endl;
    std::cout << "sum * Z_bath: " << std::setw(10) << sum_value * sr.Z_bath << std::endl;
    std::cout << "sum * Z_bath * Z_imp_correction * Z_bath_correction: " << std::setw(10)
              << sum_value * sr.Z_bath * sr.Z_imp_correction * sr.Z_bath_correction << std::endl;
  } else if (mode_name == "bare") {
    std::cout << "partition function exact: " << std::setw(10) << sr.partition_function_ref << std::endl;
    std::cout << "partition function exact * Z_imp_correction: " << std::setw(10) << sr.partition_function_ref * sr.Z_imp_correction << std::endl;
    std::cout << std::left << std::setw(10) << "order" << std::setw(30) << "value" << std::setw(30) << "time(s)" << std::setw(30)
              << "time_find_pivot(s)" << std::setw(30) << "time_pretrain(s)" << std::setw(30) << "time_train(s)" << std::endl;
    std::cout << std::setw(10) << "0" << std::setw(30) << sr.partition_function_zeroth_order_ref << std::endl;
    for (int i = 0; i < sp.order_list.size(); i++) {
      std::cout << std::setw(10) << sp.order_list[i] << std::setw(30) << sr.integral_list[i] << std::setw(30) << sr.calculation_time_list[i]
                << std::setw(30) << sr.find_pivot_time_list[i] << std::setw(30) << sr.pretrain_time_list[i] << std::setw(30) << sr.train_time_list[i]
                << std::endl;
    }
    double sum_value           = sr.partition_function_zeroth_order_ref + std::accumulate(sr.integral_list.begin(), sr.integral_list.end(), 0.0);
    double sum_time            = std::accumulate(sr.calculation_time_list.begin(), sr.calculation_time_list.end(), 0.0);
    double sum_time_find_pivot = std::accumulate(sr.find_pivot_time_list.begin(), sr.find_pivot_time_list.end(), 0.0);
    double sum_time_pretrain   = std::accumulate(sr.pretrain_time_list.begin(), sr.pretrain_time_list.end(), 0.0);
    double sum_time_train      = std::accumulate(sr.train_time_list.begin(), sr.train_time_list.end(), 0.0);
    std::cout << std::setw(10) << "sum:" << std::setw(30) << sum_value << std::setw(30) << sum_time << std::setw(30) << sum_time_find_pivot
              << std::setw(30) << sum_time_pretrain << std::setw(30) << sum_time_train << std::endl;
    std::cout << "sum * Z_imp_correction: " << std::setw(10) << sum_value * sr.Z_imp_correction << std::endl;
  } else if (mode_name == "inchworm") {
    std::cerr << "not implemented yet" << std::endl;
    std::exit(EXIT_FAILURE);
  } else {
    std::cerr << "invalid mode_name" << std::endl;
    std::exit(EXIT_FAILURE);
  }
} // end of print_summary

void ModeBase::validate_input() {
  if (gp.integral_variable == "v" && gp.tci_shape != "plain" || gp.integral_variable != "v" && gp.tci_shape == "plain") {
    std::cerr << "invalid combination of integral_variable and tci_shape" << std::endl;
    std::cerr << "integral_variable: " << gp.integral_variable << ", tci_shape: " << gp.tci_shape << std::endl;
    std::exit(EXIT_FAILURE);
  }
} // end of validate_input

void ModeBase::evaluate() {

  // setup the mapping and jacobian functions for the transformation of time-ordered variables v->tau
  cv_func change_variable;
  jb_func jacobian;
  if (tp.mapping_v == 0 || tp.mapping_v == 4) {
    change_variable = change_variable0;
    jacobian        = jacobian0;
  } else if (tp.mapping_v == 1) {
    change_variable = change_variable1;
    jacobian        = jacobian1;
  } else if (tp.mapping_v == 2) {
    change_variable = change_variable2;
    jacobian        = jacobian2;
  } else if (tp.mapping_v == 3) {
    change_variable = change_variable3;
    jacobian        = jacobian3;
  } else {
    std::cerr << "invalid mapping_v" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  for (int order : sp.order_list) {
    if (sp.debug > 0) std::cout << "order: " << order << std::endl;
    int n = 2 * order; // number of operators
    // generate valid n_left list
    std::vector<int> n_left_list(n - 1);
    std::iota(n_left_list.begin(), n_left_list.end(), 1);
    if (sp.use_bare_propagator) {
      n_left_list = {0}; // bare expansion does not need n_left
    }
    // generate valid phi and iota pairs
    std::vector<std::pair<std::vector<int>, std::vector<int>>> phi_pair_list{};
    if (!gp.do_segment) {
      std::vector<int> index_range(n);
      std::iota(index_range.begin(), index_range.end(), 0);
      phi_pair_list = get_all_phi(index_range);
    }
    // std::vector<int> index_range(n);
    // std::iota(index_range.begin(), index_range.end(), 0);
    // phi_pair_list = get_all_phi(index_range);
    std::vector<int> phi_list(phi_pair_list.size()); // phi_list is an index list, phi_pair_list contains actual phi pairs
    std::iota(phi_list.begin(), phi_list.end(), 0);
    std::vector<std::pair<std::vector<int>, std::vector<int>>> iota_pair_list{};
    if (gp.integral_variable != "v_iota") { iota_pair_list = get_all_iota(mp.gf_block_shape, order); }
    // iota_pair_list = get_all_iota(mp.gf_block_shape, order);
    std::vector<int> iota_list(iota_pair_list.size());
    std::iota(iota_list.begin(), iota_list.end(), 0); // iota_list is an index list, iota_pair_list contains actual iota pairs

    //discrete index that needed to be looped over: n_left, phi, iota (i.e., at most 3 loops); here sp.bl_index and sp.subspace_index are assumed to be fixed. In inchworm mode, these two indices are either performed with an outer loop or add to tci as an physical index
    auto loop1 = Loop("empty", std::vector<int>{0});
    auto loop2 = Loop("empty", std::vector<int>{0});
    auto loop3 = Loop("empty", std::vector<int>{0});
    if (gp.integrand == "plain" && gp.integral_variable == "v") {
      loop1 = Loop("n_left", n_left_list);
      loop2 = Loop("phi", phi_list);
      loop3 = Loop("iota", iota_list);
    } else if (gp.integrand == "plain" && gp.integral_variable == "v_iota") {
      loop1 = Loop("n_left", n_left_list);
      loop2 = Loop("phi", phi_list);
    } else if (gp.integrand == "sum_phi" && gp.integral_variable == "v") {
      loop1 = Loop("n_left", n_left_list);
      loop2 = Loop("iota", iota_list);
    } else if (gp.integrand == "sum_phi" && gp.integral_variable == "v_iota") {
      loop1 = Loop("n_left", n_left_list);
    } else {
      std::cerr << "invalid combination of integrand, integral_variable" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    auto start_time        = std::chrono::high_resolution_clock::now();
    double time_find_pivot = 0.0;
    double time_pretrain   = 0.0;
    double time_train      = 0.0;
    loop1.value            = 0;

    for (auto val1 : loop1.container) {
      loop2.value = 0;
      for (auto val2 : loop2.container) {
        loop3.value = 0;
        for (auto val3 : loop3.container) {
          int n_left = -1;
          std::vector<int> phi_d_list{};
          std::vector<int> phi_d_dag_list{};
          std::vector<int> iota_d_list{};
          std::vector<int> iota_d_dag_list{};
          int id_phi  = -1;
          int id_iota = -1;
          if (loop1.name == "n_left") {
            n_left = val1;
            std::cout << "n_left: " << n_left << std::endl;
          } else if (loop2.name == "n_left") {
            n_left = val2;
            std::cout << "n_left: " << n_left << std::endl;
          } else if (loop3.name == "n_left") {
            n_left = val3;
            std::cout << "n_left: " << n_left << std::endl;
          }
          if (loop1.name == "phi")
            id_phi = val1;
          else if (loop2.name == "phi")
            id_phi = val2;
          else if (loop3.name == "phi")
            id_phi = val3;
          if (loop1.name == "iota")
            id_iota = val1;
          else if (loop2.name == "iota")
            id_iota = val2;
          else if (loop3.name == "iota")
            id_iota = val3;
          if (id_phi != -1) {
            std::tie(phi_d_list, phi_d_dag_list) = phi_pair_list[id_phi];
          } else if (gp.integrand == "sum_phi") {
          } else {
            std::cerr << "invalid combination of integrand, integral_variable, and tci_shape" << std::endl;
            std::exit(EXIT_FAILURE);
          }
          if (id_iota != -1) {
            std::tie(iota_d_list, iota_d_dag_list) = iota_pair_list[id_iota];
          } else if (gp.integral_variable == "v_iota") {
          } else {
            std::cerr << "invalid combination of integrand, integral_variable, and tci_shape" << std::endl;
            std::exit(EXIT_FAILURE);
          }
          if (n_left == -1) {
            std::cerr << "n_left is not set" << std::endl;
            std::exit(EXIT_FAILURE);
          }

          long count     = 0;
          auto integrand = [this, &count, &phi_d_list = phi_d_list, &phi_d_dag_list = phi_d_dag_list, &iota_d_list = iota_d_list,
                            &iota_d_dag_list = iota_d_dag_list, &n_left, &change_variable, &jacobian, &phi_list,
                            &phi_pair_list](const std::vector<double> &variables) -> double {
            std::vector<double> taus_left{};
            std::vector<double> taus_right{};
            std::vector<double> taus{};
            std::vector<double> taus_left_sym{};
            std::vector<double> taus_right_sym{};
            std::vector<double> taus_sym{};
            std::vector<double> vs{};
            std::vector<double> iotas{};
            count++;
            double integrand_val = 0.0;
            // set integral variables
            if (gp.integral_variable == "v") {
              vs = variables;
            } else if (gp.integral_variable == "v_iota" && gp.tci_shape == "vertex") {
              vs.reserve(variables.size());
              iotas.reserve(variables.size());
              for (int i = 0; i < variables.size(); i++) {
                double int_part;
                double frac_part;
                frac_part = modf(variables[i], &int_part);
                vs.push_back(frac_part);
                iotas.push_back(int_part);
              }
            } else if (gp.integral_variable == "v_iota" && gp.tci_shape == "partition") {
              int mid_idx = variables.size() / 2;
              iotas       = std::vector<double>(variables.begin(), variables.begin() + mid_idx);
              vs          = std::vector<double>(variables.begin() + mid_idx, variables.end());
            } else if (gp.integral_variable == "v_iota" && gp.tci_shape == "full") {
              for (size_t i = 0; i < variables.size(); i++) {
                if (i % 2 == 0) {
                  iotas.push_back(variables[i]);
                } else {
                  vs.push_back(variables[i]);
                }
              }
            } else {
              std::cerr << "not implemented gp.integral_variable && gp.tci_shape combination" << std::endl;
              std::cerr << "gp.integral_variable: " << gp.integral_variable << ", gp.tci_shape: " << gp.tci_shape << std::endl;
              std::exit(EXIT_FAILURE);
            } // end of splitting variables
            // for bare mode, taus_left = taus_right = taus
            std::tie(taus_left, taus_right, taus) = obtain_taus(vs, n_left, sp.tau_split, sp.tau_max, change_variable);

            // if there exist duplicated element in taus, then the integrand is set to zero
            //TODO: find a more precise approximation
            if (is_duplicated(taus)) {
              std::cout << "Warning: duplicated elements in taus" << std::endl;
              if (gp.trick == "random_auxi") { return tp.auxi_height * get_hash_random_number(variables); }
              return 0.0;
            }
            // std::cout << "taus: " << std::endl;
            // print_vector(taus);
            // if (is_duplicated(taus)) {
            //   for (size_t i = 0; i < taus.size(); i++) {
            //     if (taus[i] == taus[i + 1]) {
            //       if (std::abs(taus[i] - 0) < 1e-16) {
            //         taus[i + 1] += 1e-16;
            //       } else {
            //         taus[i] -= 1e-16;
            //       }
            //     }
            //   }
            // }
            // std::cout << "taus_new: " << std::endl;
            // print_vector(taus);

            if (tp.mapping_v == 4 && sp.use_bare_propagator == false) {
              for (int i = taus_left.size() - 1; i >= 0; i--) { taus_left_sym.push_back(sp.tau_split - taus_left[i]); }
              for (int i = taus_right.size() - 1; i >= 0; i--) { taus_right_sym.push_back(sp.tau_max + sp.tau_split - taus_right[i]); }
              taus_sym = taus_left_sym;
              taus_sym.insert(taus_sym.end(), taus_right_sym.begin(), taus_right_sym.end());
            } else if (tp.mapping_v == 4 && sp.tau_split == 0.0) {
              for (int i = taus.size() - 1; i >= 0; i--) { taus_sym.push_back(sp.tau_max - taus[i]); }
            }

            std::vector<int> phi_loop_list{0};
            std::vector<std::pair<std::vector<int>, std::vector<int>>> phi_loop_pair_list{{phi_d_list, phi_d_dag_list}};
            if (gp.integrand == "sum_phi" && gp.do_segment) {
              phi_loop_pair_list = generate_phi_segment(iotas, mp.gf_block_shape);
              phi_loop_list.resize(phi_loop_pair_list.size());
              std::iota(phi_loop_list.begin(), phi_loop_list.end(), 0);
            } else if (gp.integrand == "sum_phi") {
              phi_loop_list      = phi_list;
              phi_loop_pair_list = phi_pair_list;
            }

            for (auto phi_id : phi_loop_list) {
              auto phi_d     = phi_loop_pair_list[phi_id].first;
              auto phi_d_dag = phi_loop_pair_list[phi_id].second;
              std::vector<int> iota_d{};
              std::vector<int> iota_d_dag{};
              if (gp.integral_variable == "v_iota") {
                iota_d     = get_elements_int(phi_d, iotas);
                iota_d_dag = get_elements_int(phi_d_dag, iotas);
                // sanity check
                std::vector<int> number_in_block_d     = generate_number_in_block(mp.gf_block_shape, iota_d);
                std::vector<int> number_in_block_d_dag = generate_number_in_block(mp.gf_block_shape, iota_d_dag);
                if (number_in_block_d != number_in_block_d_dag) { continue; }
              } else if (gp.integral_variable == "v") {
                iota_d     = iota_d_list;
                iota_d_dag = iota_d_dag_list;
              } else {
                std::cerr << "not implemented" << std::endl;
                std::exit(EXIT_FAILURE);
              } // end of setting iota_d and iota_d_dag
              auto tau_d         = get_elements(phi_d, taus);
              auto tau_d_dag     = get_elements(phi_d_dag, taus);
              auto integrand_phi = evaluate_u_tau_max(sr.u_tau_zeroth_order, sp.tau_split, sp.tau_max, mp.all_d_ops, mp.all_d_dag_ops,
                                                      mp.gf_block_shape, cp, mp.Delta_tau, mp.ad_imp, sr.u_interpolator, tau_d, tau_d_dag, iota_d,
                                                      iota_d_dag, sp.bl_index, sp.subspace_index, sp.use_bare_propagator, sp.gf_index, cp.gf_struct);
              double j           = jacobian(taus_right, sp.tau_max, sp.tau_split);
              if (sp.use_bare_propagator == false) { j *= jacobian(taus_left, sp.tau_split, 0.0); }
              if (tp.mapping_v == 4) {
                auto tau_d_sym     = get_elements(phi_d, taus_sym);
                auto tau_d_dag_sym = get_elements(phi_d_dag, taus_sym);
                auto integrand_phi_sym =
                   evaluate_u_tau_max(sr.u_tau_zeroth_order, sp.tau_split, sp.tau_max, mp.all_d_ops, mp.all_d_dag_ops, mp.gf_block_shape, cp,
                                      mp.Delta_tau, mp.ad_imp, sr.u_interpolator, tau_d_sym, tau_d_dag_sym, iota_d, iota_d_dag, sp.bl_index,
                                      sp.subspace_index, sp.use_bare_propagator, sp.gf_index, cp.gf_struct);
                // j_sym should be the same as j
                integrand_val += (integrand_phi_sym * j + integrand_phi * j) / 2;
              } else {
                // if (tp.mapping_v == 1 || tp.mapping_v == 3) {
                //   integrand_val += integrand_phi;
                // } else {
                //   integrand_val += integrand_phi * j;
                // }
                integrand_val += integrand_phi * j;
              }
            } // end of loop over phi
            if (gp.trick == "random_auxi") { return tp.auxi_height * get_hash_random_number(variables) + integrand_val; }
            return integrand_val;
          };

          //training
          std::vector<std::vector<double>> input{};
          std::vector<std::vector<double>> weight{};
          std::vector<int> init_pivot{};
          if (gp.integral_variable == "v") {
            for (int i = 0; i < n; i++) {
              input.push_back(tp.v_value);
              weight.push_back(tp.v_weight);
            }
          } else if (gp.integral_variable == "v_iota" && gp.tci_shape == "vertex") {
            std::vector<double> v_iota_value;
            std::vector<double> v_iota_weight;
            for (int i = 0; i < mp.n_phi; i++) {
              for (int j = 0; j < tp.v_value.size(); j++) {
                v_iota_value.push_back(i + tp.v_value[j]);
                v_iota_weight.push_back(tp.v_weight[j]);
              }
            }
            for (int i = 0; i < n; i++) {
              input.push_back(v_iota_value);
              weight.push_back(v_iota_weight);
            }
          } else if (gp.integral_variable == "v_iota" && gp.tci_shape == "partition") {
            std::vector<double> iota_value{};
            iota_value.resize(mp.n_phi);
            std::iota(iota_value.begin(), iota_value.end(), 0);
            std::vector<double> iota_weight(iota_value.size(), 1.0);
            for (int i = 0; i < n; i++) {
              input.push_back(iota_value);
              weight.push_back(iota_weight);
            }
            for (int i = 0; i < n; i++) {
              input.push_back(tp.v_value);
              weight.push_back(tp.v_weight);
            }
          } else if (gp.integral_variable == "v_iota" && gp.tci_shape == "full") {
            std::vector<double> iota_value{};
            iota_value.resize(mp.n_phi);
            std::iota(iota_value.begin(), iota_value.end(), 0);
            std::vector<double> iota_weight(iota_value.size(), 1.0);
            for (int i = 0; i < n; i++) {
              input.push_back(iota_value);
              input.push_back(tp.v_value);
              weight.push_back(iota_weight);
              weight.push_back(tp.v_weight);
            }
          } else {
            std::cerr << "not implemented" << std::endl;
            std::exit(EXIT_FAILURE);
          }

          // set initial pivot
          // smart initial pivot
          //           std::cout << "init_input: " << std::endl;
          // for (int i = 0; i < input.size(); i++) { init_pivot.push_back(int(input[i].size() / 2));}
          // std::cout << input[i][init_pivot[i]] << std::endl;
          // }

          for (int i = 0; i < input.size(); i++) { init_pivot.push_back(0); }
          // for (int i = 0; i < input.size(); i++) {
          //   for (int j = 0; j < input[i].size(); j++) {
          //     if (input[i][j] != 0) {
          //       init_pivot.push_back(j);
          //       break;
          //     }
          //     if(j == input[i].size() - 1) std::cerr << "input is all zero" << std::endl;
          //   }
          // }

          std::vector<double> init_input{};
          for (int i = 0; i < init_pivot.size(); i++) { init_input.push_back(input[i][init_pivot[i]]); }
          double init_integrand = integrand(init_input);
          if (sp.debug > 1) { std::cout << "init_integrand: " << init_integrand << std::endl; }
          if (init_integrand == 0) {
            std::cerr << "Warning: initial pivot is zero !!" << std::endl;
            std::cerr << "tau_max: " << sp.tau_max << ", tau_split: " << sp.tau_split << std::endl;
            std::cerr << "loop1.name: " << loop1.name << ", loop1.value: " << val1 << std::endl;
            std::cerr << "loop2.name: " << loop2.name << ", loop2.value: " << val2 << std::endl;
            std::cerr << "loop3.name: " << loop3.name << ", loop3.value: " << val3 << std::endl;
            continue;
          }
          std::vector<std::vector<int>> init_global_pivots{};
          if (gp.trick == "spin_pivot" && gp.integral_variable == "v_iota") {
            std::vector<int> pivot{};
            if (gp.tci_shape == "vertex") {
              for (int i = 0; i < input.size(); i++) { pivot.push_back(int(input[i].size() / 2)); }
            }
            if (gp.tci_shape == "partition") {
              for (int i = 0; i < input.size() / 2; i++) { pivot.push_back(1); }
              for (int i = input.size() / 2; i < input.size(); i++) { pivot.push_back(0); }
            }
            init_global_pivots.push_back(pivot);
          }
          // multiply jacobian for HS and AD's mapping
          double const_jacobian = 1.0;
          // if (tp.mapping_v == 1 || tp.mapping_v == 3) {
          //   std::vector<double> taus_fake_right(n - n_left, 0.0);
          //   const_jacobian = jacobian(taus_fake_right, sp.tau_max, sp.tau_split);
          //   if (sp.use_bare_propagator == false) {
          //     std::vector<double> taus_fake_left(n_left, 0.0);
          //     const_jacobian *= jacobian(taus_fake_left, sp.tau_split, 0.0);
          //   }
          // }
          // double integral = do_TCI<double, double>(integrand, input, weight, init_pivot, count, tp.sweep_bound, tp.bond_dim, tp.reltol, tp.fullPiv,
          //                                          tp.tci_prrlu, tp.error_type, tp.error_eval, tp.convergence_bound, tp.convergence_iter, sp.debug,
          //                                          init_global_pivots, const_jacobian);
          double integral = calculate_sum<double, double>(integrand, input, weight,const_jacobian);
          loop3.value += integral;
        } // end of loop3
        loop2.value += loop3.value;
      } // end of loop2
      loop1.value += loop2.value;
    } // end of loop1
    sr.integral_list.push_back(loop1.value);
    auto end_time = std::chrono::high_resolution_clock::now();
    sr.calculation_time_list.push_back(std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() / 1e6);
    sr.find_pivot_time_list.push_back(time_find_pivot);
    sr.pretrain_time_list.push_back(time_pretrain);
    sr.train_time_list.push_back(time_train);
  } // end of order loop
  std::cout << "completed" << std::endl;

} // end of evaluate_propagator

void ModeBase::evaluate_greens_function() {} // end of evaluate_green_function