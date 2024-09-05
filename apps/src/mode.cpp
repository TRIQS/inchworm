#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <h5/h5.hpp>
#include "utility.hpp"
#include "mode.hpp"

void ModeBase::clear_tci_results() {
  sr.integral_list.clear();
  sr.calculation_time_list.clear();
} // end of clear_results

void ModeBase::read_json_parameters(std::string json_file_path) {
  NVTX_RANGE("read params", 0);
  namespace pt = boost::property_tree;
  pt::ptree root;
  pt::read_json(json_file_path, root);

  // Read global parameters
  //// required parameters
  gp.target            = root.get<std::string>("gp.target");
  gp.integrand         = root.get<std::string>("gp.integrand");
  gp.integral_variable = root.get<std::string>("gp.integral_variable");
  gp.tci_shape         = root.get<std::string>("gp.tci_shape");
  gp.model_type        = root.get<int>("gp.model_type");
  //// optional parameters
  try {
    gp.ergodicity = root.get<std::string>("gp.ergodicity");
  } catch (const std::exception &e) { gp.ergodicity = "none"; }
  try {
    gp.do_segment = root.get<bool>("gp.do_segment");
  } catch (const std::exception &e) { gp.do_segment = false; }
  try {
    gp.output_prefix = root.get<std::string>("gp.output_prefix");
  } catch (const std::exception &e) { gp.output_prefix = "sim"; }
  try {
    gp.unsummed_tci = root.get<int>("gp.unsummed_tci");
  } catch (const std::exception &e) { gp.unsummed_tci = 0; }
  try {
    gp.energy_shift = root.get<double>("gp.energy_shift");
  } catch (const std::exception &e) { gp.energy_shift = 0.0; }

  // Read construction parameters
  //// required parameters
  cp.beta        = root.get<double>("cp.beta");
  cp.n_tau_green = root.get<int>("cp.n_tau_green");
  cp.n_tau_inch  = root.get<int>("cp.n_tau_inch");
  cp.n_tau_hyb   = root.get<int>("cp.n_tau_hyb");
  for (pt::ptree::value_type &g_s : root.get_child("cp.gf_struct")) {
    std::string name = g_s.first;
    int size         = g_s.second.get_value<int>();
    cp.gf_struct.emplace_back(std::make_pair(name, size));
  }

  // Read model parameters
  //// required parameters
  mp.n_site = root.get<int>("mp.n_site");
  mp.n_bath = root.get<int>("mp.n_bath");
  mp.n_spin = root.get<int>("mp.n_spin");
  mp.U      = root.get<double>("mp.U");
  mp.mu     = root.get<double>("mp.mu");
  mp.t      = root.get<double>("mp.t");
  int size  = root.get_child("mp.epsilon").size();
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
  //// required parameters
  sp.order_Chebyshev = root.get<int>("sp.order_Chebyshev");
  int debug_level    = root.get<int>("sp.debug");
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
  size = root.get_child("sp.order_list_first").size();
  sp.order_list_first.resize(size);
  i = 0;
  for (pt::ptree::value_type &order : root.get_child("sp.order_list_first")) {
    sp.order_list_first[i] = order.second.get_value<int>();
    i++;
  }
  //// optional parameters
  try {
    sp.tau_max         = root.get<double>("sp.tau_max");
    sp.tau_split_ratio = root.get<double>("sp.tau_split_ratio");
    sp.tau_split       = sp.tau_split_ratio * sp.tau_max;
    sp.bl_index        = root.get<int>("sp.bl_index");
    sp.subspace_index  = root.get<int>("sp.subspace_index");
  } catch (const std::exception &e) {
    sp.tau_max         = -1.0;
    sp.tau_split_ratio = -1.0;
    sp.tau_split       = -1.0;
    sp.bl_index        = -1.0;
    sp.subspace_index  = -1.0;
  }

  // Read TCI parameters
  //// required parameters
  tp.n_GK                 = root.get<int>("tp.n_GK");
  tp.tci_prrlu            = root.get<int>("tp.tci_prrlu");
  tp.bond_dim_init             = root.get<int>("tp.bond_dim_init");
  tp.bond_dim_increase         = root.get<int>("tp.bond_dim_increase");
  tp.bond_dim_max              = root.get<int>("tp.bond_dim_max");
  tp.sweep_bound          = root.get<int>("tp.sweep_bound");
  tp.auxi_height          = root.get<double>("tp.auxi_height");
  tp.reltol               = root.get<double>("tp.reltol");
  tp.fullPiv              = root.get<bool>("tp.fullPiv");
  tp.error_type           = root.get<int>("tp.error_type");
  tp.error_eval           = root.get<int>("tp.error_eval");
  tp.convergence_bound    = root.get<double>("tp.convergence_bound");
  tp.convergence_iter     = root.get<int>("tp.convergence_iter");
  tp.integral_lower_bound = root.get<double>("tp.integral_lower_bound");
  //// optional parameters
  try {
    tp.decay_rate = root.get<double>("tp.decay_rate");
  } catch (const std::exception &e) { tp.decay_rate = 1.0; }

  if (sp.debug > 0) { std::cout << "json parameter file read successfully" << std::endl; }
} // end of read_json_parameters

hyb_tau_t ModeBase::read_hyb_function(std::string hyb_file_path, model_params_t const &mp, constr_params_t const &cp) {

  const double TOL = 1e-14;
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
  auto Delta_tau      = hyb_tau_t{{cp.beta, Fermion, static_cast<long>(tau_grid.size())}, cp.gf_struct};
  auto Delta_tau_grid = Delta_tau[0].mesh();
  for (int i = 0; i < Delta_tau_grid.size(); i++) {
    if (std::abs(Delta_tau_grid[i] - tau_grid[i]) > TOL) {
      std::cerr << "tau_grid in hyb_file_path does not match the tau_grid generated by the code" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  for (int block = 0; block < cp.gf_struct.size(); block++) {
    for (auto [i, j] : product_range(mp.n_site, mp.n_site)) {
      std::vector<double> Delta_tau_ij{};
      Delta_tau_ij.reserve(tau_grid.size());
      h5_read(grp, "data/" + std::to_string(block) + '_' + std::to_string(i) + std::to_string(j), Delta_tau_ij);
      for (int k = 0; k < Delta_tau_grid.size(); k++) {
        Delta_tau[block][Delta_tau_grid[k]](i, j) = Delta_tau_ij[k];
        if (std::abs(Delta_tau[block](Delta_tau_grid[k])(i, j) - Delta_tau[block][Delta_tau_grid[k]](i, j)) > TOL) {
          std::cerr << "inconsistency in Delta_tau" << std::endl;
        }
      }
    }
  }

  return Delta_tau;
}

void ModeBase::prepare_input(std::string hyb_file_path) {
  NVTX_RANGE("prepare input", 0);

  // TCI setup
  std::tie(tp.v_value, tp.v_weight) = select_quadrature_GK(tp.n_GK, 0, 1);

  // operators setup
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

  // grid setup
  sp.n_tau_linear = cp.n_tau_inch;
  if (sp.order_Chebyshev != 0) {
    sp.interp_type = interpolation_type::linear_Chebyshev;
    sp.n_tot       = sp.n_tau_linear + (sp.n_tau_linear - 1) * (sp.order_Chebyshev + 1);
  } else {
    sp.interp_type = interpolation_type::cspline;
    sp.n_tot       = sp.n_tau_linear;
  }

  // model setup
  if (gp.model_type == 0) { //model_type 0: discrete bath, where exact results (reference) are available
    std::tie(sr.Z_bath_correction, sr.Z_imp_correction, sr.Z_bath, mp.Delta_tau, mp.ad_imp, sr.u_tau_ref, sr.G_tau_ref) =
       discrete_setup(mp.n_site, mp.n_bath, mp.n_spin, mp.U, mp.mu, mp.t, cp, mp.theta, mp.epsilon, sp.n_tot, sp.n_tau_linear, sp.order_Chebyshev,
                      sp.grid_linear, sp.grid);

    sr.u_interpolator_ref     = interpolator_t<scalar_t>(sr.u_tau_ref, sp.n_tot, sp.n_tau_linear, sp.order_Chebyshev, sp.interp_type, sp.grid);
    sr.partition_function_ref = trace(sr.u_interpolator_ref(cp.beta));
    if (sp.debug > 0) {
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
    mp.ad_imp                         = imp_setup(mp.n_site, mp.n_spin, mp.U, mp.mu, mp.t, cp);
    sr.Z_imp_correction               = std::exp(-(mp.ad_imp.get_gs_energy()) * cp.beta);
    mp.Delta_tau                      = read_hyb_function(hyb_file_path, mp, cp);
    std::tie(sp.grid_linear, sp.grid) = generate_inchworm_grid(0, cp.beta, sp.n_tau_linear, sp.order_Chebyshev);
    //// params that does not have analytical reference
    // sr.Z_bath = 0;
    // sr.Z_bath_correction = 0;
    // sr.u_tau_ref = {};
    // sr.G_tau_ref = {};
    // sr.partition_function_ref              = 0.0;
  } else {
    std::cerr << "invalid model_type" << std::endl;
    std::exit(EXIT_FAILURE);
  }
} // end of prepare_input

void ModeBase::print_summary() {
  if (sp.debug <= 0) return;
  // global parameters
  std::cout << "#### global params summary ####" << std::endl;
  std::cout << "mode_name: " << mode_name << std::endl;
  std::cout << "target: " << gp.target << std::endl;
  std::cout << "integrand: " << gp.integrand << std::endl;
  std::cout << "integral_variable: " << gp.integral_variable << std::endl;
  std::cout << "tci_shape: " << gp.tci_shape << std::endl;
  std::cout << "ergodicity: " << gp.ergodicity << std::endl;
  std::cout << "do_segment: " << gp.do_segment << std::endl;
  std::cout << "unsummed_tci: " << gp.unsummed_tci << std::endl;
  std::cout << "output_prefix: " << gp.output_prefix << std::endl;
  if (gp.model_type == 0) {
    std::cout << "model_type: discrete bath" << std::endl;
  } else if (gp.model_type == 1) {
    std::cout << "model_type: continuous bath" << std::endl;
  } else {
    std::cerr << "invalid model_type" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::cout << "energy_shift: " << gp.energy_shift << std::endl;
  std::cout << "#### total_time (s) ####" << std::endl;
  std::cout << sr.total_time << std::endl;
} // end of print_summary

void ModeBase::validate_input() {
  bool valid_tci_prrlu = (tp.tci_prrlu == 2);
  if (!valid_tci_prrlu) {
    std::cerr << "invalid tci_prrlu" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  bool valid_integral_variable      = (gp.integral_variable == "v_iota");
  bool valid_tci_shape              = (gp.tci_shape == "partition");
  bool valid_integrand_plus_segment = (gp.integrand == "sum_phi" && gp.do_segment) || (gp.integrand == "plain" && !gp.do_segment);

  if (!valid_integral_variable || !valid_tci_shape || !valid_integrand_plus_segment) {
    std::cerr << "invalid combination of integral_variable, tci_shape, integrand, do_segment" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  bool valid_unsumed_mode = (gp.unsummed_tci == 1 || gp.unsummed_tci == 2) && mode_name == "inchworm" || (gp.unsummed_tci == 0) && (mode_name == "bare" || mode_name == "debug");
  if (!valid_unsumed_mode) {
    std::cerr << "invalid unsummed_tci and mode_name combination" << std::endl;
    std::exit(EXIT_FAILURE);
  }
} // end of validate_input

void ModeBase::prepare_eval() {
  NVTX_RANGE("prepare eval", 0);
  // set up quantities that are valid for all calls of evaluate, i.e., those does not depend on tau and propagator
  // setup the mapping and jacobian functions for the transformation of time-ordered variables v->tau
  ep.change_variable = change_variable0;
  ep.jacobian        = jacobian0;
  if (!gp.do_segment) {
    // generate the union set elements in sp.order_list_first and sp.order_list, only store each element once
    std::vector<int> order_list_union(sp.order_list_first.size() + sp.order_list.size());
    auto it =
       std::set_union(sp.order_list_first.begin(), sp.order_list_first.end(), sp.order_list.begin(), sp.order_list.end(), order_list_union.begin());
    order_list_union.resize(it - order_list_union.begin());
    // generate phi_pair for each order
    for (int order : order_list_union) {
      std::vector<int> index_range(2 * order);
      std::iota(index_range.begin(), index_range.end(), 0);
      ep.phi_pair_order_cache[order] = get_all_phi(index_range);
    }
  }
}

void ModeBase::evaluate(std::vector<std::vector<double>> const &unsummed_input, bool is_first_interval) {
  NVTX_RANGE("evaluate", 0);
  //// set up before the order loop
  size_t unsummed_tot_size = 1;
  for (auto const &v : unsummed_input) { unsummed_tot_size *= v.size(); }

  std::vector<int> order_list = is_first_interval ? sp.order_list_first : sp.order_list;

  for (int order : order_list) {
    if (sp.debug >= 1) std::cout << "order: " << order << std::endl;
    int n_opt = 2 * order;                   // number of operators is 2*order as there are equal number of creation and annihilation operators
    std::vector<int> n_left_list(n_opt - 1); // generate valid n_left list: 1, 2, ..., n_opt-1; 0 and n_opt will not give inchworm proper diagram
    std::iota(n_left_list.begin(), n_left_list.end(), 1);
    if (sp.use_bare_propagator) {
      n_left_list = {0}; // bare expansion does not need n_left
    }

    //For the case of no segment: according to order, read the corresponding phi_pair_list from phi_pair_order_cache
    const std::vector<std::pair<std::vector<int>, std::vector<int>>> *phi_pair_list = nullptr;
    std::vector<int> phi_list = {}; // phi_list is an index list, phi_pair_list holds actual phi pairs
    if (!gp.do_segment) {
      auto it = ep.phi_pair_order_cache.find(order);
      if (it == ep.phi_pair_order_cache.end()) {
        std::cerr << "phi_pair_order_cache does not have the order: " << order << std::endl;
        std::exit(EXIT_FAILURE);
      }
      phi_pair_list = &(it->second);
      phi_list      = std::vector<int>(phi_pair_list->size());
      std::iota(phi_list.begin(), phi_list.end(), 0);
    }

    //discrete index that needed to be looped over: n_left, phi; iota is assumed to be in the tensor train as physical indices
    auto loop1 = Loop("empty", std::vector<int>{0});
    auto loop2 = Loop("empty", std::vector<int>{0});
    if (gp.integrand == "plain" && gp.integral_variable == "v_iota") {
      loop1 = Loop("n_left", n_left_list);
      loop2 = Loop("phi", phi_list);
    } else if (gp.integrand == "sum_phi" && gp.integral_variable == "v_iota") {
      loop1 = Loop("n_left", n_left_list);
    } else {
      std::cerr << "invalid combination of integrand, integral_variable" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // loop started
    auto start_time = std::chrono::high_resolution_clock::now();
    loop1.value     = std::vector<double>(unsummed_tot_size, 0);
    for (auto val1 : loop1.container) {
      loop2.value = std::vector<double>(unsummed_tot_size, 0);
      for (auto val2 : loop2.container) {

        // set up quantities
        double auxi_height = tp.auxi_height;
        int n_left         = -1;
        std::vector<int> phi_d_list{};
        std::vector<int> phi_d_dag_list{};
        int id_phi = -1;
        if (loop1.name == "n_left") {
          n_left = val1;
        } else if (loop2.name == "n_left") {
          n_left = val2;
        }
        if (sp.debug > 1) { std::cout << "n_left: " << n_left << std::endl; }
        if (loop1.name == "phi")
          id_phi = val1;
        else if (loop2.name == "phi")
          id_phi = val2;
        if (id_phi != -1) {
          std::tie(phi_d_list, phi_d_dag_list) = (*phi_pair_list)[id_phi];
        } else if (gp.integrand == "sum_phi") {
        } else {
          std::cerr << "invalid combination of integrand, integral_variable, and tci_shape" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (gp.integral_variable == "v_iota") {
        } else {
          std::cerr << "invalid combination of integrand, integral_variable, and tci_shape" << std::endl;
          std::exit(EXIT_FAILURE);
        }

        long count     = 0;
        auto integrand = [this, &count, &phi_d_list = phi_d_list, &phi_d_dag_list = phi_d_dag_list,  &n_left, &auxi_height](std::vector<double> variables) -> double {
          NVTX_RANGE("integrand", 0);
          std::vector<double> taus_left{};
          std::vector<double> taus_right{};
          std::vector<double> taus{};
          std::vector<double> vs{};
          std::vector<double> iotas{};
          count++;

          int bl_index       = sp.bl_index;
          int subspace_index = sp.subspace_index;
          double tau_max     = sp.tau_max;
          if ((gp.unsummed_tci == 1 || gp.unsummed_tci == 2)) {
            int index = static_cast<int>(variables[0]);
            std::tie(bl_index, subspace_index) = bl1_to_bl2(index, mp.ad_imp.get_subspace_dims());
          }

          if (gp.unsummed_tci == 2) {
            tau_max = variables[1];
          }

          // set integral variables
          if (gp.integral_variable == "v_iota" && gp.tci_shape == "partition") {
            int mid_idx = variables.size() / 2;
            iotas       = std::vector<double>(variables.begin() + gp.unsummed_tci, variables.begin()+ gp.unsummed_tci + mid_idx);
            vs          = std::vector<double>(variables.begin() + gp.unsummed_tci + mid_idx, variables.end());
          } else {
            std::cerr << "invalid combination of integral_variable and tci_shape" << std::endl;
            std::exit(EXIT_FAILURE);
          }

          // for bare mode, taus_left = taus_right = taus
          std::tie(taus_left, taus_right, taus) = obtain_taus(vs, n_left, sp.tau_split, tau_max, ep.change_variable);
          if (is_duplicated(taus)) {
            std::cerr << "Warning: duplicated elements in taus" << std::endl;
            if (gp.ergodicity == "random_auxi_adaptive") { return auxi_height * get_random_number(); }
            return 0.0;
          }
          if (if_contains(taus, sp.tau_split)) {
            std::cerr << "Warning: tau_split is in taus" << std::endl;
            if (gp.ergodicity == "random_auxi_adaptive") { return auxi_height * get_random_number(); }
            return 0.0;
          }

          std::vector<int> phi_loop_list{0};
          std::vector<std::pair<std::vector<int>, std::vector<int>>> phi_loop_pair_list{{phi_d_list, phi_d_dag_list}};
          if (gp.integrand == "sum_phi" && gp.do_segment) {
            NVTX_RANGE("segment", 5);
            auto it = ep.phi_pair_iota_cache.find(iotas);
            if (it == ep.phi_pair_iota_cache.end()) {
              phi_loop_pair_list            = generate_phi_segment(iotas, mp.gf_block_shape);
              ep.phi_pair_iota_cache[iotas] = phi_loop_pair_list;
            } else {
              phi_loop_pair_list = it->second;
            }
            phi_loop_list.resize(phi_loop_pair_list.size());
            std::iota(phi_loop_list.begin(), phi_loop_list.end(), 0);
          } 
          
          double integrand_val = 0.0;
          for (auto phi_id : phi_loop_list) { // possible phi loop due to sum_phi
            auto phi_d     = phi_loop_pair_list[phi_id].first;
            auto phi_d_dag = phi_loop_pair_list[phi_id].second;
            std::vector<int> iota_d{};
            std::vector<int> iota_d_dag{};
            if (gp.integral_variable == "v_iota") {
              iota_d     = get_elements_int(phi_d, iotas);
              iota_d_dag = get_elements_int(phi_d_dag, iotas);
              std::vector<int> number_in_block_d     = generate_number_in_block(mp.gf_block_shape, iota_d);
              std::vector<int> number_in_block_d_dag = generate_number_in_block(mp.gf_block_shape, iota_d_dag);
              if (number_in_block_d != number_in_block_d_dag) { continue; } // this does not give a non-zero contribution
            } else {
              std::cerr << "only v_iota is implemented" << std::endl;
              std::exit(EXIT_FAILURE);
            } // end of setting iota_d and iota_d_dag
            auto tau_d         = get_elements(phi_d, taus);
            auto tau_d_dag     = get_elements(phi_d_dag, taus);
            auto integrand_phi = evaluate_diagram(sr.u_tau_zeroth_order, sp.tau_split, tau_max, mp.all_d_ops, mp.all_d_dag_ops, mp.gf_block_shape,
                                                    cp, mp.Delta_tau, mp.ad_imp, sr.u_interpolator, tau_d, tau_d_dag, iota_d, iota_d_dag, bl_index,
                                                    subspace_index, sp.use_bare_propagator, sp.gf_index, cp.gf_struct, gp.energy_shift);
            double j           = ep.jacobian(taus_right, tau_max, sp.tau_split);
            if (sp.use_bare_propagator == false) { j *= ep.jacobian(taus_left, sp.tau_split, 0.0); }
            integrand_val += integrand_phi * j;
          } // end of loop over phi

          if (gp.ergodicity == "random_auxi_adaptive") {
            // std::cout << auxi_height * get_random_number() + integrand_val << std::endl;
            return auxi_height * get_random_number() + integrand_val;
          }
          // std::cout << integrand_val << std::endl;
          return integrand_val;
        }; // end of integrand lambda function

        //set up input, weight
        std::vector<std::vector<double>> input{};
        std::vector<std::vector<double>> weight{};
        std::vector<int> init_pivot{};
        if (gp.integral_variable == "v_iota" && gp.tci_shape == "partition") {
          std::vector<double> iota_value{};
          iota_value.resize(mp.n_phi);
          std::iota(iota_value.begin(), iota_value.end(), 0);
          std::vector<double> iota_weight(iota_value.size(), 1.0);
          for (int i = 0; i < n_opt; i++) {
            input.push_back(iota_value);
            weight.push_back(iota_weight);
          }
          for (int i = 0; i < n_opt; i++) {
            input.push_back(tp.v_value);
            weight.push_back(tp.v_weight);
          }
        } else {
          std::cerr << "not implemented" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (gp.unsummed_tci != 0 && mode_name == "inchworm") {
          for (int i = unsummed_input.size() - 1; i >= 0; i--) {
            input.insert(input.begin(), unsummed_input[i]);
            weight.insert(weight.begin(), std::vector<double>(unsummed_input[i].size(), 1.0));
          }
        }

        // set up initial pivot, initial input, initial integrand
        for (int i = 0; i < input.size(); i++) { init_pivot.push_back(0); }
        std::vector<double> init_input{};
        for (int i = 0; i < init_pivot.size(); i++) { init_input.push_back(input[i][init_pivot[i]]); }
        double init_integrand = integrand(init_input);
        if (sp.debug > 1) { std::cout << "init_integrand: " << init_integrand << std::endl; }
        if (init_integrand == 0) {
          std::cerr << "Warning: initial integrand is zero !!" << std::endl;
          std::cerr << "tau_max: " << sp.tau_max << ", tau_split: " << sp.tau_split << std::endl;
          std::cerr << "loop1.name: " << loop1.name << ", loop1.name->val: " << val1 << std::endl;
          std::cerr << "loop2.name: " << loop2.name << ", loop2.name->val: " << val2 << std::endl;
          std::cerr << std::endl;
          continue;
        }

        
        std::vector<std::vector<int>> init_global_pivots{};
        double const_jacobian = 1.0;
        std::vector<double> integral(unsummed_tot_size, 0.0);

        bool adaptive_error = false;
        if (gp.ergodicity == "random_auxi_adaptive") { adaptive_error = true; }
        integral =
           do_TCI<double, double>(integrand, input, weight, init_pivot, count, tp.sweep_bound, tp.bond_dim_init, tp.bond_dim_increase, tp.bond_dim_max, tp.reltol, tp.fullPiv, tp.tci_prrlu,
                                  tp.error_type, tp.error_eval, tp.convergence_bound, tp.convergence_iter, sp.debug, init_global_pivots,
                                  const_jacobian, gp.unsummed_tci, unsummed_tot_size, adaptive_error, tp.decay_rate, &auxi_height);

        for (size_t i = 0; i < unsummed_tot_size; i++) { loop1.value[i] += integral[i]; }
      } // end of loop2
    } // end of loop1
    sr.integral_list.push_back(loop1.value);
    auto end_time = std::chrono::high_resolution_clock::now();
    sr.calculation_time_list.push_back(std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() / 1e6);
  } // end of order loop
  sr.total_time += std::accumulate(sr.calculation_time_list.begin(), sr.calculation_time_list.end(), 0.0);

} // end of evaluate
