#include "./insert.hpp"
#include "../diagram/print.hpp"

namespace inchworm::moves {

  mc_weight_t insert::attempt() {
    // Pick up the value of alpha and choose the operators
    int n_fops = (data->h_diag.get_fops()).size();
    int li1    = rng(n_fops);
    int li2    = rng(n_fops);

    // Choice of times for insertion.
    double tau1 = rng(data->tau_max());
    double tau2 = rng(data->tau_max());

    int N        = data->size();
    bool success = false;
    while (not success) success = data->insert(tau1, li1, tau2, li2);

    if (data->use_bare_propagator) {
      auto diagram = data->get_time_diagram();
      print_configuration(diagram);
      new_w_hyb        = determinant(diagram);
      auto new_U_frame = propagator_product(data->h_diag, diagram, data->tau_max());
      new_w_loc        = new_U_frame.frobenius_norm();
    } else {
      std::vector<double> split_times = {data->tau_split()};
      auto diagram                    = data->get_time_diagram(split_times);
      new_w_hyb                       = inclusion_exclusion(diagram);
      auto new_U_frame                = propagator_product(data->U_tau, data->h_diag, diagram, data->tau_max());

      new_w_loc = new_U_frame.frobenius_norm();
    }
    
    if(std::abs(new_w_hyb)<1e-10) exit(0);

    auto w_hyb_ratio = new_w_hyb / data->last_accepted_w_hyb;

    // atomic weight
    auto w_loc_ratio = new_w_loc / data->last_accepted_w_loc;

    // proposition probability
    auto t_ratio = std::pow(data->tau_max() / (N + 1), 2);

    //printf("salut0:  %f  \n", t_ratio * w_loc_ratio * w_hyb_ratio);
    if (verbose > 1) {
      printf("new_w_hyb, new_w_loc, last_w_hyb, last_w_loc, t_ratio:  %f %f   %f %f   %f\n", new_w_hyb, new_w_loc, data->last_accepted_w_hyb, data->last_accepted_w_loc, t_ratio);
      fflush(stdout);
    }
    return t_ratio * w_loc_ratio * w_hyb_ratio;
  }

  mc_weight_t insert::accept() {

    data->last_accepted_w_hyb = new_w_hyb;
    data->last_accepted_w_loc = new_w_loc;
    //data->last_accepted_U_frame = new_U_frame;

    return 1;
  }

  void insert::reject() { data->erase_last(); }

} // namespace inchworm::moves
