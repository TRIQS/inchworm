#include "./qmc_config.hpp"

namespace inchworm {

  qmc_config_t::qmc_config_t(params_t const &params, atom_diag const &h_diag, block_gf_const_view<imtime> delta)
     : tau_max_(params.beta), h_diag(h_diag), delta(delta), last_accepted_U_frame(h_diag) {

    auto propagator_struct = find_propagator_struct(h_diag);
    U_tau                  = u_tau_t{{params.beta, Fermion, params.n_tau}, propagator_struct};
    assign_identity_to_propagator(U_tau, 0); // assign identity matrices to the frame 0 of U_tau

    last_accepted_w_loc = 1.0;
    last_accepted_w_hyb = 1.0;
  }

  double qmc_config_t::tau_max() { return tau_max_; }
  double qmc_config_t::tau_split() { return tau_split_; }
  int qmc_config_t::size() { return c_list.size(); }

  bool qmc_config_t::try_insert(double tau, int linear_index, double tau_dag, int linear_index_dag) {
    update_lists();
    for (int i = 0; i < size() - 1; i++)
      if ((c_list[i].tau == tau) or (cdag_list[i].tau == tau_dag)) return false;
    c_list.push_back({tau, linear_index});
    cdag_list.push_back({tau_dag, linear_index_dag});
    return true;
  }
  bool qmc_config_t::try_erase(int i, int i_dag) {
    update_lists();
    if ((size() <= i) or (size() <= i_dag)) return false;
    c_list.erase(c_list.begin() + i);
    cdag_list.erase(cdag_list.begin() + i_dag);
    return true;
  }
  void qmc_config_t::update_lists() {
    c_list = last_accepted_c_list;
    cdag_list = last_accepted_cdag_list;
  }
  void qmc_config_t::clear() {
    c_list.clear();
    cdag_list.clear();
  }

  time_diagram_t qmc_config_t::get_time_diagram(std::vector<double> const &split_times) { return time_diagram_t(c_list, cdag_list, split_times); }

  time_diagram_t qmc_config_t::get_time_diagram() {
    std::vector<double> split_times{};
    return time_diagram_t(c_list, cdag_list, split_times);
  }
} // namespace inchworm
