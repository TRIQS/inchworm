#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "./mode.hpp"

using namespace inchworm;

void ModeNestedTCI::read_json_parameters(std::string json_file_path) {

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

  // Read TCI parameters for iota
  tp_iota.tci_prrlu            = root.get<bool>("tp_iota.tci_prrlu");
  tp_iota.bond_dim             = root.get<int>("tp_iota.bond_dim");
  tp_iota.sweep_bound          = root.get<int>("tp_iota.sweep_bound");
  tp_iota.integral_error_bound = root.get<double>("tp_iota.integral_error_bound");
  tp_iota.pivot_error_bound    = root.get<double>("tp_iota.pivot_error_bound");
}

void ModeNestedTCI::init(std::string json_file_path) { ModeNestedTCI::read_json_parameters(json_file_path); }

void ModeNestedTCI::run_single_element() {

  std::vector<int> iotai(mp.n_phi);
  std::iota(iotai.begin(), iotai.end(), 0);
  auto wi_iota = std::vector(mp.n_phi, 1.0);
  for (int order : sp.order_list) {
    auto start_time = std::chrono::high_resolution_clock::now();
    int n           = 2 * order; // number of tau's
    std::vector<int> pivot1(n, 0);
    std::vector<int> index_range(n);
    std::iota(index_range.begin(), index_range.end(), 0);
    auto phi_pair_list = get_all_phi(index_range); //gives all possible phi

    std::vector<int> iota_pivots(n, 0);           // this is an intermediate variable for generating all possible iota
    std::vector<int> iota_pivots_range(mp.n_phi); // the int version of iotai
    std::iota(iota_pivots_range.begin(), iota_pivots_range.end(), 0);
    std::vector<std::vector<int>> all_iota_pivots{};
    generate_combinations(iota_pivots_range, iota_pivots, 0, all_iota_pivots);

    double integral_sum_phi = 0.0;
    for (auto [phi_d_list, phi_d_dag_list] : phi_pair_list) {
      double integral_sum_n_left = 0.0;
      for (int n_left = 1; n_left < n; n_left++) {
        double integral_sum_iota        = 0.0;
        long count_iota                 = 0;
        auto get_u_tau_max_element_iota = [this, &count_iota, &phi_d_list = phi_d_list, &phi_d_dag_list = phi_d_dag_list, &n_left, &pivot1,
                                           &n](const std::vector<int> &iotas) {
          std::vector<int> iota_d_list     = get_elements(phi_d_list, iotas);
          std::vector<int> iota_d_dag_list = get_elements(phi_d_dag_list, iotas);
          long count                       = 0;
          auto get_u_tau_max_element       = [this, &count, &phi_d_list = phi_d_list, &phi_d_dag_list = phi_d_dag_list, &iota_d_list = iota_d_list,
                                        &iota_d_dag_list = iota_d_dag_list, &n_left](const std::vector<double> &vs) {
            auto [taus_left, taus_right, taus] = obtain_taus(vs, n_left, sp.tau_split, sp.tau_max);
            double integrand = evaluate_u_tau_max(sr.u_tau_max_zeroth_order, sp.tau_split, sp.tau_max, mp.all_d_ops, mp.all_d_dag_ops,
                                                        mp.gf_block_shape, cp, mp.Delta_tau, mp.ad_imp, sr.u_interpolator, get_elements(phi_d_list, taus),
                                                        get_elements(phi_d_dag_list, taus), iota_d_list, iota_d_dag_list, sp.bl_index, sp.subspace_index);
            count++;
            double j = jacobian(taus_left, sp.tau_split, 0.0) * jacobian(taus_right, sp.tau_max, sp.tau_split);
            return integrand * j;
          };
          std::vector<double> vs1;
          for (int i = 0; i < pivot1.size(); i++) { vs1.push_back(tp.vi[pivot1[i]]); }
          double u_tau_max_element_vs1 = get_u_tau_max_element(vs1);
          if (sp.debug > 1) {
            auto [taus_left1, taus_right1, taus1] = obtain_taus(vs1, n_left, sp.tau_split, sp.tau_max);
            print_pivot1(iota_d_list, iota_d_dag_list, get_elements(phi_d_list, taus1), get_elements(phi_d_dag_list, taus1), u_tau_max_element_vs1);
          }
          if (u_tau_max_element_vs1 == 0) {
            count_iota++;
            return 0.0;
          }
          auto input              = std::vector(n, tp.vi);
          auto weight             = std::vector(n, tp.wi_v);
          double integral_element = do_TCI<double, double>(get_u_tau_max_element, input, weight, pivot1, tp.sweep_bound, tp.bond_dim,
                                                           tp.integral_error_bound, tp.pivot_error_bound, tp.tci_prrlu, sp.debug, count);
          return integral_element;
        }; // end of definition of get_u_tau_max_element_iota

        std::vector<int> pivot1_iota = {};
        for (auto iota_pivot1 : all_iota_pivots) {
          auto value = get_u_tau_max_element_iota(iota_pivot1);
          if (value != 0) {
            pivot1_iota = iota_pivot1;
            break;
          }
        }
        if (pivot1_iota.empty()) { continue; }
        auto input_iota         = std::vector(n, iotai);
        auto weight_iota        = std::vector(n, wi_iota);
        double integral_element = do_TCI<double, int>(get_u_tau_max_element_iota, input_iota, weight_iota, pivot1_iota, tp_iota.sweep_bound, tp_iota.bond_dim,
                                                         tp_iota.integral_error_bound, tp_iota.pivot_error_bound, tp_iota.tci_prrlu, sp.debug, count_iota);
        integral_sum_iota += integral_element;
        integral_sum_n_left += integral_sum_iota;
      }
      integral_sum_phi += integral_sum_n_left;
    }
    auto end_time            = std::chrono::high_resolution_clock::now();
    auto duration            = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    auto duration_in_seconds = static_cast<double>(duration) / 1e6;
    sr.calculation_time_list.push_back(duration_in_seconds);
    sr.integral_order_list.push_back(integral_sum_phi);
  }
}
