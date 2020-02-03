#include "./qmc_data.hpp"

namespace inchworm {

  bool config_t::try_insert(double tau, int linear_index, double tau_dag, int linear_index_dag) {
    //std::printf("try_erase? %d   %2.4f %d  %2.4f %d \n", size(),  tau, linear_index, tau_dag, linear_index_dag);
    for (int i = 0; i < size() - 1; i++)
      if ((c_list[i].tau == tau) or (cdag_list[i].tau == tau_dag)) return false;
    c_list.push_back({tau, linear_index}); // FIXME upper_bound ordering
    cdag_list.push_back({tau_dag, linear_index_dag});
    return true;
  }

  bool config_t::try_erase(int i, int i_dag) {
    //std::printf("try_erase? %d   %d %d \n", size(), i, i_dag);
    if ((size() <= i) or (size() <= i_dag)) return false;
    c_list.erase(c_list.begin() + i);
    cdag_list.erase(cdag_list.begin() + i_dag);
    return true;
  }


  qmc_data_t::qmc_data_t(params_t const &params, atom_diag const &h_diag, u_tau_t const &u_tau, block_gf_const_view<imtime> delta)
     : w{1.0, 1.0}, tau_max(params.beta), h_diag(h_diag), delta(delta), u_frame(h_diag), u_tau(u_tau) {

    //auto propagator_struct = find_propagator_struct(h_diag);
    //u_tau                  = u_tau_t{{params.beta, Fermion, params.n_tau}, propagator_struct};
    //assign_identity_to_propagator(u_tau, 0); // assign identity matrices to the frame 0 of u_tau
  }

} // namespace inchworm
