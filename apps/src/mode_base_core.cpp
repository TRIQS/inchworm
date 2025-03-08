#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <h5/h5.hpp>
#include "utility.hpp"
#include "mode.hpp"
#include "save.hpp"

void ModeBase::prepare_input() {
  NVTX_RANGE("prepare input", 0);

  // TCI setup
  if (gp.map_type == 0) {
    std::tie(tp.v_value, tp.v_weight) = select_quadrature_GK(tp.n_GK, 0, 1);
  } else if (gp.map_type == 3) {
    int n_grid_tanh_sinh  = tp.n_GK;
    double h_tanh_sinh    = 4.0 / n_grid_tanh_sinh;
    auto [xi_tanh_sinh, wi_tanh_sinh] = tanh_sinh_quadrature(0, 1, h_tanh_sinh, n_grid_tanh_sinh);
    double tol            = 1e-14;
    tp.v_value  = xi_tanh_sinh;
    tp.v_weight = wi_tanh_sinh;
  } else {
    std::cerr << "Invalid map_type" << std::endl;
    std::exit(EXIT_FAILURE);
  }

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
    if (sp.debug > 0 && rank == 0) {
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
    mp.Delta_tau                      = read_hyb_function();
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

void ModeBase::validate_input() {
  //check the consistency among bl1_to_bl3, bl1_to_bl2, bl2_to_bl1, bl2_to_bl3, bl3_to_bl1
  std::vector<int> block_shape = mp.ad_imp.get_subspace_dims();
  int bl_index                 = int(block_shape.size() / 2);
  int subspace_index           = block_shape[bl_index] - 1;
  int i_index                  = 0;
  int j_index                  = block_shape[bl_index] - 1;
  int total_index              = 0;
  for (int i = 0; i < bl_index; i++) { total_index += block_shape[i] * block_shape[i]; }
  total_index += i_index * block_shape[bl_index] + j_index;

  // bl1_to_bl3
  auto [bl, i, j] = bl1_to_bl3(total_index, block_shape);
  if (bl != bl_index || i != i_index || j != j_index) {
    std::cerr << "bl: " << bl << ", i: " << i << ", j: " << j << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // bl1_to_bl2
  auto [bl_x, subspace_x] = bl1_to_bl2(total_index, block_shape);
  if (bl_x != bl_index || subspace_x != subspace_index) {
    std::cerr << "bl: " << bl_x << ", subspace: " << subspace_x << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // bl2_to_bl1
  int total_index_ = bl2_to_bl1(bl_index, subspace_index, block_shape);
  if (total_index_ != total_index) {
    std::cout << "total_index: " << total_index << std::endl;
    std::cerr << "total_index_: " << total_index_ << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // bl2_to_bl3
  auto [bl_, i_, j_] = bl2_to_bl3(bl_index, subspace_index, block_shape);
  if (bl_ != bl_index || i_ != i_index || j_ != j_index) {
    std::cerr << "bl_: " << bl_ << ", i_: " << i_ << ", j_: " << j_ << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // bl3_to_bl1
  auto [bl_index___, total_index__] = bl3_to_bl1(bl_index, i_index, j_index, block_shape);
  if (total_index__ != total_index) {
    std::cerr << "total_index__: " << total_index__ << std::endl;
    std::exit(EXIT_FAILURE);
  }

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

  bool valid_unsumed_mode_propagator = (gp.unsummed_tci == 1 || gp.unsummed_tci == 2) && mode_name == "inchworm"
     || (gp.unsummed_tci == 0) && (mode_name == "bare" || mode_name == "debug");
  bool valid_unsumed_mode_gf = mode_name == "inchworm" || mode_name == "debug";
  if (!(valid_unsumed_mode_propagator && gp.target == "propagator") && !(valid_unsumed_mode_gf && gp.target == "greens_function") && !(valid_unsumed_mode_gf && gp.target == "greens_function_restart") ) {
    std::cerr << "invalid unsummed_tci and mode_name combination" << std::endl;
    std::exit(EXIT_FAILURE);
  }
} // end of validate_input

void ModeBase::prepare_eval() {
  NVTX_RANGE("prepare eval", 0);
  // set up quantities that are valid for all calls of evaluate, i.e., those does not depend on tau and propagator
  // setup the mapping and jacobian functions for the transformation of time-ordered variables v->tau
  //set 0
  if (gp.map_type == 0) {
    ep.change_variable = change_variable0;
    ep.jacobian        = jacobian0;
  } else if (gp.map_type == 3) {
    // // set 3
    ep.change_variable = change_variable3;
    ep.jacobian        = jacobian3;
  } else {
    std::cerr << "Invalid map_type" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (gp.target == "greens_function" && gp.do_segment) {
    throw std::runtime_error("gp.target == greens_function and gp.do_segment == true is not supported yet");
  }
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
  // generate the space for sr.statistics
  if (mode_name == "inchworm") {
    sr.statistics.resize(sp.n_tau_linear - 1);
  } else {
    sr.statistics.resize(1);
  }
}

void ModeBase::evaluate(std::vector<std::vector<double>> const &unsummed_input, bool is_first_interval, size_t inchworm_index) {
  NVTX_RANGE("evaluate", 0);
  //// set up before the order loop
  size_t unsummed_tot_size = 1;
  for (auto const &v : unsummed_input) { unsummed_tot_size *= v.size(); }

  std::vector<int> order_list = is_first_interval ? sp.order_list_first : sp.order_list;

  std::vector<std::tuple<int, int, int>> looptotl;
  std::vector<Loop<std::vector<int>>> loop1_list;
  std::vector<Loop<std::vector<int>>> loop2_list;
  for (size_t order_ind = 0; order_ind < order_list.size(); order_ind++) {
    int order = order_list[order_ind];
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
    for (size_t loop1_idx = 0; loop1_idx < loop1.container.size(); loop1_idx++) {
      for (size_t loop2_idx = 0; loop2_idx < loop2.container.size(); loop2_idx++) {
        looptotl.push_back(std::make_tuple(order_ind, loop1_idx, loop2_idx));
      }
    }
    loop1_list.push_back(loop1);
    loop2_list.push_back(loop2);
  }

  std::vector<double> integral_order_flatten(order_list.size() * unsummed_tot_size, 0);
  std::vector<double> integral_order_rank_flatten(order_list.size() * unsummed_tot_size, 0);

  /// statistics start
  std::vector<long> func_evals_order(order_list.size(), 0);
  std::vector<long> func_evals_order_rank(order_list.size(), 0);
  std::vector<long> warning_same_time_order(order_list.size(), 0);
  std::vector<long> warning_same_time_order_rank(order_list.size(), 0);
  std::vector<long> warning_tau_split_order(order_list.size(), 0);
  std::vector<long> warning_tau_split_order_rank(order_list.size(), 0);
  std::vector<long> warning_tau_max_order(order_list.size(), 0);
  std::vector<long> warning_tau_max_order_rank(order_list.size(), 0);
  std::vector<double> max_diff_order(order_list.size(), 0);
  std::vector<double> max_diff_order_rank(order_list.size(), 0);
  std::vector<double> max_auxi_height_order(order_list.size(), 0);
  std::vector<double> max_auxi_height_order_rank(order_list.size(), 0);
  std::vector<double> max_error_order(order_list.size(), 0);
  std::vector<double> max_error_order_rank(order_list.size(), 0);
  std::vector<double> u_tau_sum_order(order_list.size(), 0);
  std::vector<double> u_tau_sum_order_rank(order_list.size(), 0);
  std::vector<double> time_order(order_list.size(), 0);
  std::vector<double> time_order_rank(order_list.size(), 0);
  std::vector<long> nTCI_order(order_list.size(), 0);
  std::vector<long> nTCI_order_rank(order_list.size(), 0);
  std::vector<double> integral_max_order(order_list.size(), 0);
  std::vector<double> integral_max_order_rank(order_list.size(), 0);
  /// statistics end

  // evenly distribute the work among ranks, save the index for current rank in a vector current_loop_index, thoses index will be assigned to the current_loop_index one by one
  size_t num_loops = looptotl.size();
  for (size_t ii = 0; ii < num_loops; ii++) {
    auto index      = looptotl[ii];
    auto start_time = std::chrono::high_resolution_clock::now();
    int order_idx   = std::get<0>(index);
    int loop1_idx   = std::get<1>(index);
    int loop2_idx   = std::get<2>(index);
    if (ii % size != rank) { continue; }
    if (sp.debug >= 3) {
      std::cout << "rank: " << rank << ", order_idx: " << order_idx << ", loop1_idx: " << loop1_idx << ", loop2_idx: " << loop2_idx << std::endl;
    }
    int order  = order_list[order_idx];
    int n_opt  = 2 * order;
    auto loop1 = loop1_list[order_idx];
    auto loop2 = loop2_list[order_idx];
    auto val1  = loop1.container[loop1_idx];
    auto val2  = loop2.container[loop2_idx];

    //For the case of no segment: according to order, read the corresponding phi_pair_list from phi_pair_order_cache
    const std::vector<std::pair<std::vector<int>, std::vector<int>>> *phi_pair_list = nullptr;
    if (!gp.do_segment) {
      auto it = ep.phi_pair_order_cache.find(order);
      if (it == ep.phi_pair_order_cache.end()) {
        std::cerr << "phi_pair_order_cache does not have the order: " << order << std::endl;
        std::exit(EXIT_FAILURE);
      }
      phi_pair_list = &(it->second);
    }

    //       // set up quantities
    double auxi_height = tp.auxi_height;
    int seed           = 0;
    int n_left         = -1;
    std::vector<int> phi_d_list{};
    std::vector<int> phi_d_dag_list{};
    int id_phi = -1;
    if (loop1.name == "n_left") {
      n_left = val1;
    } else if (loop2.name == "n_left") {
      n_left = val2;
    }
    if (sp.debug > 1 && rank == 0) { std::cout << "n_left: " << n_left << std::endl; }
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

    std::string debug_info;
    {
      std::ostringstream oss;
      oss << "rank: " << rank << ", order: " << order << ", n_left: " << n_left << ", tau_split: " << std::setprecision(16) << sp.tau_split;
      debug_info = oss.str();
    }

    long count         = 0;
    double mini_height = 0.0;
    double mini_value  = 1e10;
    auto integrand     = [this, &count, &phi_d_list = phi_d_list, &phi_d_dag_list = phi_d_dag_list, &n_left, &auxi_height, &seed, &mini_height,
                      &mini_value, &debug_info, &order_idx, &warning_same_time_order_rank, &warning_tau_split_order_rank,
                      &warning_tau_max_order_rank](std::vector<double> variables) -> double {
      NVTX_RANGE("integrand", 0);
      std::vector<double> taus_left{};
      std::vector<double> taus_right{};
      std::vector<double> taus{};
      std::vector<double> vs{};
      std::vector<double> iotas{};
      count++;

      int bl_index              = sp.bl_index;
      int subspace_index        = sp.subspace_index;
      double tau_max            = sp.tau_max;
      double tau_split          = sp.tau_split;
      std::vector<int> gf_index = sp.gf_index;
      if (gp.target == "propagator") {
        if ((gp.unsummed_tci == 1 || gp.unsummed_tci == 2)) {
          int index                          = static_cast<int>(variables[0]);
          std::tie(bl_index, subspace_index) = bl1_to_bl2(index, mp.ad_imp.get_subspace_dims());
        }
        if (gp.unsummed_tci == 2) { tau_max = variables[1]; }
      } else if (gp.target == "greens_function") {
        if (gp.unsummed_tci == 1) {
          tau_split = variables[0];
        } else if (gp.unsummed_tci == 2) {
          int index       = int(variables[0]);
          tau_split       = variables[1];
          auto [bl, i, j] = bl1_to_bl3(index, mp.gf_block_shape);
          gf_index[0]     = bl2_to_bl1_gf(bl, i, mp.gf_block_shape);
          gf_index[1]     = bl2_to_bl1_gf(bl, j, mp.gf_block_shape);
        }
      }

      // set integral variables
      if (gp.integral_variable == "v_iota" && gp.tci_shape == "partition") {
        int mid_idx = (variables.size() - gp.unsummed_tci) / 2;
        iotas       = std::vector<double>(variables.begin() + gp.unsummed_tci, variables.begin() + gp.unsummed_tci + mid_idx);
        vs          = std::vector<double>(variables.begin() + gp.unsummed_tci + mid_idx, variables.end());
      } else {
        std::cerr << "invalid combination of integral_variable and tci_shape" << std::endl;
        std::exit(EXIT_FAILURE);
      }

      // for bare mode, taus_left = taus_right = taus
      std::tie(taus_left, taus_right, taus) = obtain_taus(vs, n_left, tau_split, tau_max, ep.change_variable);
      // adjust taus
      std::tie(taus_left, taus_right, taus) = obtain_taus_restricted(taus_left, taus_right, taus, 1e-15, n_left, tau_split, tau_max);

      if (is_duplicated(taus)) {
        // mini_height = 0.0;
        std::cerr << "Warning: duplicated elements in taus" << std::endl;
        std::cerr << "taus: ";
        for (auto tau : taus) { std::cerr << tau << " "; }
        std::cerr << std::endl;
        std::cerr << "debug_info: " << debug_info << std::endl;
        warning_same_time_order_rank[order_idx] += 1;
        //gracifally exit the program
        // std::exit(EXIT_FAILURE);
        // if (gp.ergodicity == "random_auxi_adaptive") { return auxi_height * get_random_number() + mini_height; }
        // if (gp.ergodicity == "random_auxi_adaptive") { return auxi_height * get_hash_random_number(variables) + mini_height; }
        if (gp.ergodicity == "random_auxi_adaptive") { return 0.0; }
        return 0.0;
      }

      if (if_contains(taus, tau_split)) {
        // mini_height = 0.0;
        std::cerr << "Warning: tau_split is in taus" << std::endl;
        std::cerr << "taus: ";
        for (auto tau : taus) { std::cerr << tau << " "; }
        std::cerr << std::endl;
        std::cerr << "tau_split: " << tau_split << std::endl;
        std::cerr << "debug_info: " << debug_info << std::endl;
        warning_tau_split_order_rank[order_idx] += 1;
        // std::exit(EXIT_FAILURE);
        // if (gp.ergodicity == "random_auxi_adaptive") { return auxi_height * get_random_number() + mini_height; }
        // if (gp.ergodicity == "random_auxi_adaptive") { return auxi_height * get_hash_random_number(variables) + mini_height; }
        if (gp.ergodicity == "random_auxi_adaptive") { return 0.0; }
        return 0.0;
      }

      if (if_contains(taus, tau_max)) {
        // mini_height = 0.0;
        std::cerr << "Warning: tau_max is in taus" << std::endl;
        std::cerr << "taus: ";
        for (auto tau : taus) { std::cerr << tau << " "; }
        std::cerr << std::endl;
        std::cerr << "tau_max: " << tau_max << std::endl;
        std::cerr << "debug_info: " << debug_info << std::endl;
        warning_tau_max_order_rank[order_idx] += 1;
        // std::exit(EXIT_FAILURE);
        // if (gp.ergodicity == "random_auxi_adaptive") { return auxi_height * get_random_number() + mini_height; }
        // if (gp.ergodicity == "random_auxi_adaptive") { return auxi_height * get_hash_random_number(variables)+ mini_height; }
        if (gp.ergodicity == "random_auxi_adaptive") { return 0.0; }
        return 0.0;
      }

      std::vector<int> phi_loop_list{0};
      std::vector<std::pair<std::vector<int>, std::vector<int>>> phi_loop_pair_list{{phi_d_list, phi_d_dag_list}};
      if (gp.integrand == "sum_phi" && gp.do_segment) {
        NVTX_RANGE("segment", 5);
        if (gp.do_cache) {
          auto it = ep.phi_pair_iota_cache.find(iotas);
          if (it == ep.phi_pair_iota_cache.end()) {
            phi_loop_pair_list            = generate_phi_segment(iotas, mp.gf_block_shape);
            ep.phi_pair_iota_cache[iotas] = phi_loop_pair_list;
          } else {
            phi_loop_pair_list = it->second;
          }
        } else {
          phi_loop_pair_list = generate_phi_segment(iotas, mp.gf_block_shape);
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
          iota_d                                 = get_elements_int(phi_d, iotas);
          iota_d_dag                             = get_elements_int(phi_d_dag, iotas);
          std::vector<int> number_in_block_d     = generate_number_in_block(mp.gf_block_shape, iota_d);
          std::vector<int> number_in_block_d_dag = generate_number_in_block(mp.gf_block_shape, iota_d_dag);
          if (number_in_block_d != number_in_block_d_dag) { continue; } // this does not give a non-zero contribution
        } else {
          std::cerr << "only v_iota is implemented" << std::endl;
          std::exit(EXIT_FAILURE);
        } // end of setting iota_d and iota_d_dag
        auto tau_d     = get_elements(phi_d, taus);
        auto tau_d_dag = get_elements(phi_d_dag, taus);
        auto integrand_phi =
           evaluate_diagram(sr.u_tau_zeroth_order, tau_split, tau_max, mp.all_d_ops, mp.all_d_dag_ops, mp.gf_block_shape, cp, mp.Delta_tau, mp.ad_imp,
                                sr.u_interpolator, tau_d, tau_d_dag, iota_d, iota_d_dag, bl_index, subspace_index, sp.use_bare_propagator, gf_index,
                                cp.gf_struct, mp.theta, mp.epsilon, mp.n_bath, gp.model_type, gp.exponent_u, gp.do_enum);
        double j = ep.jacobian(taus_right, tau_max, tau_split);
        if (sp.use_bare_propagator == false) { j *= ep.jacobian(taus_left, tau_split, 0.0); }
        integrand_val += integrand_phi * j;
      } // end of loop over phi

      if (abs(integrand_val) < mini_value && abs(integrand_val) > 1e-16) { mini_value = integrand_val; }

      // if (gp.ergodicity == "random_auxi_adaptive") { return auxi_height * get_random_number() + integrand_val; }
      // if (gp.ergodicity == "random_auxi_adaptive") { return auxi_height * get_hash_random_number(variables) + integrand_val; }
      if (gp.ergodicity == "random_auxi_adaptive") { return integrand_val; }
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

      std::vector<double> v_value_left{};
      std::vector<double> v_value_right{};
      std::vector<double> v_weight_left{};
      std::vector<double> v_weight_right{};
      if (gp.do_adaptive_nGK) {
        int nGK_left = tp.n_GK;
        if (n_left > 0) { int nGK_left = adjust_nGK(n_left, 3, tp.n_GK, sp.tau_split, 0.0, ep.change_variable); }
        int nGK_right = adjust_nGK(n_opt - n_left, 3, tp.n_GK, sp.tau_max, sp.tau_split, ep.change_variable);
        if (sp.debug > 0) {
          std::cout << "-----adaptive nGK-----" << std::endl;
          std::cout << "tau_max: " << sp.tau_max << std::endl;
          std::cout << "tau_split: " << sp.tau_split << std::endl;
          std::cout << "n_opt: " << n_opt << std::endl;
          std::cout << "n_left: " << n_left << std::endl;
          std::cout << "nGK_left: " << nGK_left << std::endl;
          std::cout << "nGK_right: " << nGK_right << std::endl;
        }
        std::tie(v_value_left, v_weight_left)   = select_quadrature_GK(nGK_left, 0, 1);
        std::tie(v_value_right, v_weight_right) = select_quadrature_GK(nGK_right, 0, 1);
      } else {
        v_value_left   = tp.v_value;
        v_value_right  = tp.v_value;
        v_weight_left  = tp.v_weight;
        v_weight_right = tp.v_weight;
      }

      for (int i = 0; i < n_left; i++) {
        input.push_back(v_value_left);
        weight.push_back(v_weight_left);
      }
      for (int i = n_left; i < n_opt; i++) {
        input.push_back(v_value_right);
        weight.push_back(v_weight_right);
      }

    } else {
      std::cerr << "not implemented" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (gp.unsummed_tci != 0) {
      for (int i = unsummed_input.size() - 1; i >= 0; i--) {
        input.insert(input.begin(), unsummed_input[i]);
        weight.insert(weight.begin(), std::vector<double>(unsummed_input[i].size(), 1.0));
      }
    }

    // std::cout << "mini_value (old)" << mini_value << std::endl;
    // set up initial pivot, initial input, initial integrand
    for (int i = 0; i < input.size(); i++) { init_pivot.push_back(0); }
    std::vector<double> init_input{};
    for (int i = 0; i < init_pivot.size(); i++) { init_input.push_back(input[i][init_pivot[i]]); }
    double init_integrand = integrand(init_input);

    if (init_integrand == 0 && gp.ergodicity != "random_auxi_adaptive") {
      // try all possible pivots for the first element
      for (int i = 0; i < input[0].size(); i++) {
        init_pivot[0] = i;
        init_input[0] = input[0][init_pivot[0]];
        init_integrand = integrand(init_input);
        if (init_integrand != 0) { 
          std::cerr << "update init_pivot[0] to " << i << std::endl;
          break; }
      }
      if (init_integrand == 0) {
        std::cerr << "Warning: initial integrand is zero !!" << std::endl;
        std::cerr << "debug_info: " << debug_info << std::endl;
        continue;
      }
    }

    if (sp.debug > 1 && rank == 0) { std::cout << "init_integrand: " << init_integrand << std::endl; }

    std::vector<std::vector<int>> init_global_pivots{};
    if (gp.do_global_pivot) {
      std::vector<int> global_pivots_last_element{};
      for (int i = 0; i < input.size(); i++) { global_pivots_last_element.push_back(input[i].size() - 1); }
      init_global_pivots.push_back(global_pivots_last_element);
    }
    double const_jacobian = 1.0;
    std::vector<double> integral(unsummed_tot_size, 0.0);
    bool adaptive_error = false;
    if (gp.ergodicity == "random_auxi_adaptive") { adaptive_error = true; }
    integral = do_TCI<double, double>(
       rank, debug_info, integrand, input, weight, init_pivot, count, tp.sweep_bound, tp.bond_dim_init, tp.bond_dim_increase, tp.bond_dim_max,
       tp.reltol, tp.fullPiv, tp.tci_prrlu, tp.error_type, tp.error_eval, tp.convergence_bound, tp.convergence_iter, sp.debug, init_global_pivots,
       const_jacobian, gp.unsummed_tci, unsummed_tot_size, adaptive_error, tp.decay_rate, max_diff_order_rank, max_auxi_height_order_rank,
       max_error_order_rank, nTCI_order_rank, integral_max_order_rank, order_idx, &auxi_height, &seed, &mini_height, &mini_value);
    // std::cout << "integral" << integral[0] << std::endl;
    // std::cout << "mini_value" << mini_value << std::endl;
    // accumulate statistics
    func_evals_order_rank[order_idx] += count;
    // warning_same_time, warning_tau_split, warning_tau_max are already accumulated in the integrand lambda function
    // max_diff, max_auxi_height, max_pivot_error are already accumulated in do_TCI
    for (size_t i = 0; i < unsummed_tot_size; i++) { u_tau_sum_order_rank[order_idx] += std::abs(integral[i]); }
    auto end_time = std::chrono::high_resolution_clock::now();
    time_order_rank[order_idx] += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() / 1e6;
    // nTCI, integral_max are already accumulated in do_TCI

    // accumate integral
    for (size_t i = 0; i < unsummed_tot_size; i++) { integral_order_rank_flatten[order_idx * unsummed_tot_size + i] += integral[i]; }
  } // end of loop over looptotl

  // gather the integrals from all ranks
  MPI_Allreduce(integral_order_rank_flatten.data(), integral_order_flatten.data(), integral_order_flatten.size(), MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);

  // gather statistics
  MPI_Allreduce(func_evals_order_rank.data(), func_evals_order.data(), func_evals_order_rank.size(), MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(warning_same_time_order_rank.data(), warning_same_time_order.data(), warning_same_time_order_rank.size(), MPI_LONG, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(warning_tau_split_order_rank.data(), warning_tau_split_order.data(), warning_tau_split_order_rank.size(), MPI_LONG, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(warning_tau_max_order_rank.data(), warning_tau_max_order.data(), warning_tau_max_order_rank.size(), MPI_LONG, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(max_diff_order_rank.data(), max_diff_order.data(), max_diff_order_rank.size(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(max_auxi_height_order_rank.data(), max_auxi_height_order.data(), max_auxi_height_order_rank.size(), MPI_DOUBLE, MPI_MAX,
                MPI_COMM_WORLD);
  MPI_Allreduce(max_error_order_rank.data(), max_error_order.data(), max_error_order_rank.size(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(u_tau_sum_order_rank.data(), u_tau_sum_order.data(), u_tau_sum_order_rank.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(time_order_rank.data(), time_order.data(), time_order_rank.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(nTCI_order_rank.data(), nTCI_order.data(), nTCI_order_rank.size(), MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(integral_max_order_rank.data(), integral_max_order.data(), integral_max_order_rank.size(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  // transform flat integral to nested integral, each element is a vector of fixed order
  for (size_t i = 0; i < order_list.size(); i++) {
    sr.integral_list.push_back(
       std::vector<double>(integral_order_flatten.begin() + i * unsummed_tot_size, integral_order_flatten.begin() + (i + 1) * unsummed_tot_size));
  }
  // store statistics
  sr.statistics[inchworm_index].func_evals_order        = func_evals_order;
  sr.statistics[inchworm_index].warning_same_time_order = warning_same_time_order;
  sr.statistics[inchworm_index].warning_tau_split_order = warning_tau_split_order;
  sr.statistics[inchworm_index].warning_tau_max_order   = warning_tau_max_order;
  sr.statistics[inchworm_index].max_diff_order          = max_diff_order;
  sr.statistics[inchworm_index].max_auxi_height_order   = max_auxi_height_order;
  sr.statistics[inchworm_index].max_error_order         = max_error_order;
  sr.statistics[inchworm_index].u_tau_sum_order         = u_tau_sum_order;
  sr.statistics[inchworm_index].time_order              = time_order;
  sr.statistics[inchworm_index].nTCI_order              = nTCI_order;
  sr.statistics[inchworm_index].integral_max_order      = integral_max_order;
} // end of evaluate

std::tuple<u_tau_t, double> ModeBase::regularize_propagator(const u_tau_t &u, const model_params_t &mp, const simulation_params_t &sp,
                                                            size_t idx_tau_max, double amplification) {
  auto u_regularized = u;
  double u_max       = -1.0;
  // find the maximum element in u at idx_tau_max
  for (auto bl : range(mp.ad_imp.n_subspaces())) {
    for (auto i : range(mp.ad_imp.get_subspace_dim(bl))) {
      for (auto j : range(mp.ad_imp.get_subspace_dim(bl))) { u_max = std::max(u_max, std::abs(u[bl][idx_tau_max](i, j))); }
    }
  }
  // calculate the exponent
  double exponent = std::log(amplification * u_max) / sp.grid[idx_tau_max];
  for (auto bl : range(mp.ad_imp.n_subspaces())) {
    for (auto i : range(mp.ad_imp.get_subspace_dim(bl))) {
      for (auto j : range(mp.ad_imp.get_subspace_dim(bl))) {
        for (size_t i_tau = 0; i_tau < sp.grid.size(); i_tau++) { u_regularized[bl][i_tau](i, j) *= std::exp(-exponent * sp.grid[i_tau]); }
      }
    }
  }

  return {u_regularized, exponent};
}

void ModeBase::evaluate_greens_function_bold() {
  gp.unsummed_tci = gp.unsummed_tci_green_function; // Fixme later, not elegant
  if (gp.do_regularization) {
    double exponent              = 0.0;
    std::tie(sr.u_tau, exponent) = regularize_propagator(sr.u_tau, mp, sp, sp.grid.size() - 1, gp.amplification_u);
    gp.exponent_u += exponent;
  }
  sr.u_interpolator      = interpolator_t<scalar_t>(sr.u_tau, sp.n_tot, sp.n_tau_linear, sp.order_Chebyshev, sp.interp_type, sp.grid);
  sp.use_bare_propagator = false;
  // Calculate Tr U(beta)
  scalar_t Tr_Ubeta = 0.0;
  for (int bl = 0; bl < sr.u_tau.size(); bl++) Tr_Ubeta += trace(sr.u_tau[bl][sp.n_tot - 1]);
  std::cout << "Z = Tr[U(beta)]: " << Tr_Ubeta << std::endl;

  sr.G_tau     = g_tau_t{{cp.beta, Fermion, cp.n_tau_green}, cp.gf_struct};
  int tot_dims = 0;
  for (auto subspace_dim : mp.gf_block_shape) { tot_dims += subspace_dim * subspace_dim; }
  // treat the first and last point separately
  EXPECTS(sr.G_tau[0].mesh()[0] == 0.0);
  frame_t g_frame_0 = make_bare_g_frame(mp.ad_imp, sr.u_tau, cp.gf_struct, 0.0, cp.beta) / Tr_Ubeta;
  set_frame(g_frame_0, sr.G_tau, 0);
  EXPECTS(sr.G_tau[0].mesh()[cp.n_tau_green - 1] == cp.beta);
  frame_t g_frame_beta = make_bare_g_frame(mp.ad_imp, sr.u_tau, cp.gf_struct, cp.beta, cp.beta) / Tr_Ubeta;
  set_frame(g_frame_beta, sr.G_tau, cp.n_tau_green - 1);

  //set zero-th order
  for (size_t n = 1; n < cp.n_tau_green - 1; n++) {
    sp.tau_split            = sr.G_tau[0].mesh()[n];
    auto frame_zeroth_order = make_bare_g_frame(mp.ad_imp, sr.u_tau, cp.gf_struct, sp.tau_split, cp.beta);
    for (int bl = 0; bl < mp.gf_block_shape.size(); bl++) {
      for (auto orb_d : range(mp.gf_block_shape[bl])) {
        for (auto orb_ddag : range(mp.gf_block_shape[bl])) {
          int orb_d_index                  = bl2_to_bl1_gf(bl, orb_d, mp.gf_block_shape);
          int orb_ddag_index               = bl2_to_bl1_gf(bl, orb_ddag, mp.gf_block_shape);
          sr.G_tau[bl][n](orb_d, orb_ddag) = (frame_zeroth_order[bl](orb_d, orb_ddag)) / Tr_Ubeta;
        }
      }
    }
  }
  if (gp.unsummed_tci == 0) {
    sp.tau_max = cp.beta;
    for (size_t n = 1; n < cp.n_tau_green - 1; n++) {
      std::cout << "evaluating tau[" << n << "] = " << sr.G_tau[0].mesh()[n] << std::endl;
      sp.tau_split = sr.G_tau[0].mesh()[n];
      // the most naive way is to loop over the spin-oribital index for d and ddag separately, However, since the Green's function is saved in block format, we only allow spin-orbital indices within the same block.
      for (int bl = 0; bl < mp.gf_block_shape.size(); bl++) {
        for (auto orb_d : range(mp.gf_block_shape[bl])) {
          for (auto orb_ddag : range(mp.gf_block_shape[bl])) {
            int orb_d_index    = bl2_to_bl1_gf(bl, orb_d, mp.gf_block_shape);
            int orb_ddag_index = bl2_to_bl1_gf(bl, orb_ddag, mp.gf_block_shape);
            sp.gf_index.clear();
            sp.gf_index.push_back(orb_d_index);
            sp.gf_index.push_back(orb_ddag_index);
            ModeBase::clear_tci_results();
            ModeBase::evaluate();
            double total_integral = 0.0;
            for (auto integral : sr.integral_list) total_integral += std::accumulate(integral.begin(), integral.end(), 0.0);
            double tci_result = total_integral / Tr_Ubeta;
            sr.G_tau[bl][n](orb_d, orb_ddag) += tci_result;
          }
        }
      }
    }
  } else if (gp.unsummed_tci == 1) {
    // include the time index into the tci
    std::vector<std::vector<double>> unsummed_input;
    std::vector<double> tau_list;
    for (size_t n = 1; n < cp.n_tau_green - 1; n++) { tau_list.push_back(sr.G_tau[0].mesh()[n]); }
    unsummed_input.push_back(tau_list);
    size_t dims_tau = tau_list.size();
    sp.tau_max      = cp.beta;
    sp.tau_split    = sr.G_tau[0].mesh()[1];
    for (int bl = 0; bl < mp.gf_block_shape.size(); bl++) {
      for (auto orb_d : range(mp.gf_block_shape[bl])) {
        for (auto orb_ddag : range(mp.gf_block_shape[bl])) {
          int orb_d_index    = bl2_to_bl1_gf(bl, orb_d, mp.gf_block_shape);
          int orb_ddag_index = bl2_to_bl1_gf(bl, orb_ddag, mp.gf_block_shape);
          sp.gf_index.clear();
          sp.gf_index.push_back(orb_d_index);
          sp.gf_index.push_back(orb_ddag_index);
          ModeBase::clear_tci_results();
          ModeBase::evaluate(unsummed_input);
          // save results
          std::vector<double> integrals(dims_tau, 0.0);
          for (auto integral_order : sr.integral_list) {
            for (size_t i = 0; i < integral_order.size(); i++) { integrals[i] += integral_order[i]; }
          }
          for (size_t n = 1; n < cp.n_tau_green - 1; n++) {
            auto tci_result = integrals[n - 1] / Tr_Ubeta;
            sr.G_tau[bl][n](orb_d, orb_ddag) += tci_result;
          }
        }
      }
    }
  } else if (gp.unsummed_tci == 2) {
    sp.gf_index.clear();
    sp.gf_index.push_back(0);
    sp.gf_index.push_back(0);
    //include the two orb indexs and tau index into tensor train
    std::vector<std::vector<double>> unsummed_input;
    size_t dims_orb = 0;
    for (int bl = 0; bl < mp.gf_block_shape.size(); bl++) { dims_orb += mp.gf_block_shape[bl] * mp.gf_block_shape[bl]; }
    std::vector<double> input(dims_orb, 0.0);
    std::iota(input.begin(), input.end(), 0);
    unsummed_input.push_back(input);
    std::vector<double> tau_list;
    for (size_t n = 1; n < cp.n_tau_green - 1; n++) { tau_list.push_back(sr.G_tau[0].mesh()[n]); }
    size_t dims_tau = tau_list.size();
    unsummed_input.push_back(tau_list);
    size_t dims_tot = dims_orb * dims_tau;
    sp.tau_max      = cp.beta;
    ModeBase::clear_tci_results();
    ModeBase::evaluate(unsummed_input);
    std::vector<double> integrals(dims_tot, 0.0);
    for (auto integral_order : sr.integral_list) {
      for (size_t i = 0; i < integral_order.size(); i++) { integrals[i] += integral_order[i]; }
    }
    for (size_t tau_index = 0; tau_index < dims_tau; tau_index++) {
      for (size_t orb_index = 0; orb_index < dims_orb; orb_index++) {
        size_t k                = tau_index + dims_tau * orb_index;
        auto [bl, i, j]         = bl1_to_bl3(orb_index, mp.gf_block_shape);
        size_t tau_index_actual = tau_index + 1;
        auto tci_result         = integrals[k] / Tr_Ubeta;
        sr.G_tau[bl][tau_index_actual](i, j) += tci_result;
      }
    }
  } else {
    throw std::runtime_error("evaluate Green's function: invalid unsummed_tci");
  }

  if (rank == 0) {
    // print the result
    for (size_t n = 0; n < cp.n_tau_green; n++) {
      for (int bl = 0; bl < mp.gf_block_shape.size(); bl++) {
        for (auto orb_d : range(mp.gf_block_shape[bl])) {
          for (auto orb_ddag : range(mp.gf_block_shape[bl])) {
            std::cout << "G[" << n << "][" << bl << "][" << orb_d << "][" << orb_ddag << "] = " << sr.G_tau[bl][n](orb_d, orb_ddag) << std::endl;
            if(gp.model_type==0){
            std::cout << "G_ref[" << n << "][" << bl << "][" << orb_d << "][" << orb_ddag << "] = " << sr.G_tau_ref[bl][n](orb_d, orb_ddag)
                      << std::endl;}
          }
        }
      }
    }
    auto file_name = gp.output_prefix + ".h5";
    h5::file file{file_name, 'a'};
    h5::group group{file};
    h5_save_params(this, group, "params");
    h5_save_gf(this, group, "gf", sr.G_tau);
    if(gp.model_type==0){
    h5_save_gf(this, group, "gf_ref", sr.G_tau_ref);}
  }
}
