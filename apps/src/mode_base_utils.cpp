#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <h5/h5.hpp>
#include "utility.hpp"
#include "mode.hpp"
#include "save.hpp"

void ModeBase::clear_tci_results() { sr.integral_list.clear(); } // end of clear_results

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
    gp.do_regularization = root.get<bool>("gp.do_regularization");
  } catch (const std::exception &e) { gp.do_regularization = true; }
  try {
    gp.amplification_u = root.get<double>("gp.amplification_u");
    if (gp.amplification_u < 0) { throw std::invalid_argument("Error: amplification_u should be non-negative"); }
  } catch (const std::exception &e) { gp.amplification_u = 1.0; }
  try {
    gp.do_cache = root.get<bool>("gp.do_cache");
  } catch (const std::exception &e) { gp.do_cache = true; }
  try {
    gp.do_adaptive_nGK = root.get<bool>("gp.do_adaptive_nGK");
  } catch (const std::exception &e) { gp.do_adaptive_nGK = false; }
  try {
    gp.do_global_pivot = root.get<bool>("gp.do_global_pivot");
  } catch (const std::exception &e) { gp.do_global_pivot = false; }
  try {
    gp.map_type = root.get<int>("gp.map_type");
  } catch (const std::exception &e) { gp.map_type = 0; }
  try {
    gp.do_enum = root.get<bool>("gp.do_enum");
  } catch (const std::exception &e) { gp.do_enum = false; }

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
  try {
    mp.rescale = root.get<double>("mp.rescale");
  } catch (const std::exception &e) { mp.rescale = 1.0; }
  // rescale beta, U, mu, t, epsilon, theta
  cp.beta *= mp.rescale;
  mp.U /= mp.rescale;
  mp.mu /= mp.rescale;
  mp.t /= mp.rescale;
  for (auto &e : mp.epsilon) { e /= mp.rescale; }
  for (auto &e : mp.theta) { e /= mp.rescale; }

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
    sp.inch_start_index = root.get<int>("sp.inch_start_index");
    sp.inch_end_index   = root.get<int>("sp.inch_end_index");
    sp.tau_max          = root.get<double>("sp.tau_max");
    sp.tau_split_ratio  = root.get<double>("sp.tau_split_ratio");
    sp.tau_split        = sp.tau_split_ratio * sp.tau_max;
    sp.bl_index         = root.get<int>("sp.bl_index");
    sp.subspace_index   = root.get<int>("sp.subspace_index");
  } catch (const std::exception &e) {
    sp.inch_start_index = 1;
    sp.inch_end_index   = sp.n_tau_linear - 1;
    sp.tau_max          = -1.0;
    sp.tau_split_ratio  = -1.0;
    sp.tau_split        = -1.0;
    sp.bl_index         = -1.0;
    sp.subspace_index   = -1.0;
  }

  // Read TCI parameters
  //// required parameters
  tp.n_GK              = root.get<int>("tp.n_GK");
  tp.tci_prrlu         = root.get<int>("tp.tci_prrlu");
  tp.bond_dim_init     = root.get<int>("tp.bond_dim_init");
  tp.bond_dim_increase = root.get<int>("tp.bond_dim_increase");
  tp.bond_dim_max      = root.get<int>("tp.bond_dim_max");
  tp.sweep_bound       = root.get<int>("tp.sweep_bound");
  tp.auxi_height       = root.get<double>("tp.auxi_height");
  tp.reltol            = root.get<double>("tp.reltol");
  tp.fullPiv           = root.get<bool>("tp.fullPiv");
  tp.error_type        = root.get<int>("tp.error_type");
  tp.error_eval        = root.get<int>("tp.error_eval");
  tp.convergence_bound = root.get<double>("tp.convergence_bound");
  tp.convergence_iter  = root.get<int>("tp.convergence_iter");
  //// optional parameters
  try {
    tp.decay_rate = root.get<double>("tp.decay_rate");
  } catch (const std::exception &e) { tp.decay_rate = 1.0; }

  if (sp.debug > 0 && rank == 0) { std::cout << "json parameter file read successfully" << std::endl; }
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

void ModeBase::print_params() {
  if (sp.debug <= 2) return;
  std::cout << "#### global params ####" << std::endl;
  std::cout << "target: " << gp.target << std::endl;
  std::cout << "integrand: " << gp.integrand << std::endl;
  std::cout << "integral_variable: " << gp.integral_variable << std::endl;
  std::cout << "tci_shape: " << gp.tci_shape << std::endl;
  std::cout << "ergodicity: " << gp.ergodicity << std::endl;
  std::cout << "model_type: " << gp.model_type << std::endl;
  std::cout << "do_segment: " << gp.do_segment << std::endl;
  std::cout << "do_cache: " << gp.do_cache << std::endl;
  std::cout << "output_prefix: " << gp.output_prefix << std::endl;
  std::cout << "unsummed_tci: " << gp.unsummed_tci << std::endl;
  std::cout << "do_regularization: " << gp.do_regularization << std::endl;
  std::cout << "amplification_u: " << gp.amplification_u << std::endl;
  std::cout << "exponent_u: " << gp.exponent_u << std::endl;
  std::cout << "Z_energy_shift_correction: " << gp.Z_energy_shift_correction << std::endl;
  std::cout << "do_adaptive_nGK: " << gp.do_adaptive_nGK << std::endl;
  std::cout << "do_global_pivot: " << gp.do_global_pivot << std::endl;
  std::cout << "map_type: " << gp.map_type << std::endl;
  std::cout << "do_enum: " << gp.do_enum << std::endl;

  std::cout << "#### model params ####" << std::endl;
  std::cout << "n_site: " << mp.n_site << std::endl;
  std::cout << "n_bath: " << mp.n_bath << std::endl;
  std::cout << "n_spin: " << mp.n_spin << std::endl;
  std::cout << "n_phi: " << mp.n_phi << std::endl;
  std::cout << "rescale: " << mp.rescale << std::endl;
  std::cout << "U: " << mp.U << std::endl;
  std::cout << "mu: " << mp.mu << std::endl;
  std::cout << "t: " << mp.t << std::endl;

  std::cout << "#### tci params ####" << std::endl;
  std::cout << "n_GK: " << tp.n_GK << std::endl;
  std::cout << "tci_prrlu: " << tp.tci_prrlu << std::endl;
  std::cout << "bond_dim_init: " << tp.bond_dim_init << std::endl;
  std::cout << "bond_dim_increase: " << tp.bond_dim_increase << std::endl;
  std::cout << "bond_dim_max: " << tp.bond_dim_max << std::endl;
  std::cout << "sweep_bound: " << tp.sweep_bound << std::endl;
  std::cout << "auxi_height: " << tp.auxi_height << std::endl;
  std::cout << "reltol: " << tp.reltol << std::endl;
  std::cout << "fullPiv: " << tp.fullPiv << std::endl;
  std::cout << "error_type: " << tp.error_type << std::endl;
  std::cout << "error_eval: " << tp.error_eval << std::endl;
  std::cout << "convergence_bound: " << tp.convergence_bound << std::endl;
  std::cout << "convergence_iter: " << tp.convergence_iter << std::endl;
  std::cout << "decay_rate: " << tp.decay_rate << std::endl;

  std::cout << "#### simulation params ####" << std::endl;
  std::cout << "debug: " << sp.debug << std::endl;
  std::cout << "inch_start_index: " << sp.inch_start_index << std::endl;
  std::cout << "inch_end_index: " << sp.inch_end_index << std::endl;
  std::cout << "tau_max: " << sp.tau_max << std::endl;
  std::cout << "tau_split: " << sp.tau_split << std::endl;
  std::cout << "tau_split_ratio: " << sp.tau_split_ratio << std::endl;
  std::cout << "n_tau_linear: " << sp.n_tau_linear << std::endl;
  std::cout << "order_Chebyshev: " << sp.order_Chebyshev << std::endl;
  std::cout << "n_tot: " << sp.n_tot << std::endl;
  std::cout << "bl_index: " << sp.bl_index << std::endl;
  std::cout << "subspace_index: " << sp.subspace_index << std::endl;
  std::cout << "use_bare_propagator: " << sp.use_bare_propagator << std::endl;
  std::cout << "eval_type: " << sp.eval_type << std::endl;
  std::cout << "n_skip: " << sp.n_skip << std::endl;

} // end of print_params