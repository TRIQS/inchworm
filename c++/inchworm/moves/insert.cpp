#include "./insert.hpp"
#include "../diagram/print.hpp"

namespace inchworm::moves {

  mc_weight_t insert::attempt() {

    // temporary values:
    new_config = data.config; // we first copy last accepted config before proposing the new insert
    new_w      = data.w;

    //print_diag(data.get_time_diagram());
    int n_fops = (data.h_diag.get_fops()).size();
    int li     = rng(n_fops);
    int li_dag = rng(n_fops);

    double tau     = rng(data.tau_max);
    double tau_dag = rng(data.tau_max);
    if (not new_config.try_insert(tau, li, tau_dag, li_dag)) return 0;

    //std::printf("success? %d \n", success);
    //}

    //std::printf("\ninserting:");
    if (data.use_bare_propagator) {
      auto diagram = diagram::time_diagram_t{new_config.c_list, new_config.cdag_list, {}}; // make a free function (not member of data)
      //print_configuration(diagram);
      new_w._hyb       = 0; // diagram::determinant(diagram);
      auto new_u_frame = propagator_product(data.h_diag, diagram, data.tau_max);
      new_w_loc        = new_u_frame.frobenius_norm();
    } else {
      auto diagram = diagram::time_diagram_t{new_config.c_list, new_config.cdag_list, {data.tau_split}}; // make a free function (not member of data)
      new_w._hyb   = 0;                                                                                  //diagram::inclusion_exclusion(diagram);
      auto new_u_frame = propagator_product(data.u_tau, data.h_diag, diagram, data.tau_max);
      new_w_loc        = new_u_frame.frobenius_norm();
    }

    auto w_hyb_ratio = new_w._hyb / data.w._hyb;
    auto w_loc_ratio = new_w._loc / data.w._loc;
    auto t_ratio     = std::pow(data.tau_max / (data.size() + 1), 2);

    //printf("helloaaaaaaaaaaiii02:  %f  \n", t_ratio * w_loc_ratio * w_hyb_ratio);
    if (false) {
      printf("new_w_hyb, new_w_loc, last_w_hyb, last_w_loc, t_ratio:  %f %f   %f %f   %f\n", new_w._hyb, new_w._loc, data.w._hyb, data.w._loc,
             t_ratio);
      fflush(stdout);
    }
    return t_ratio * w_loc_ratio * w_hyb_ratio;
  }

  //
  mc_weight_t insert::accept() {
    //std::printf("yes\n");
    //print_configuration(data.get_time_diagram()  );
    data.w._hyb = new_w._hyb;
    data.w._loc = new_w._loc;

    return 1;
  }

} // namespace inchworm::moves
