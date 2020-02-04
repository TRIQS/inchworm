#include "./insert.hpp"
#include "../diagram/print.hpp"

namespace inchworm::moves {

  scalar_t insert::attempt() {

    proposed_config = data.config; // we first copy last accepted config before proposing the new insert
    proposed_w      = data.w;

    int n_fops = (data.h_diag.get_fops()).size();
    int li     = rng(n_fops);
    int li_dag = rng(n_fops);

    double tau     = rng(data.tau_max);
    double tau_dag = rng(data.tau_max);
    if (not proposed_config.try_insert(tau, li, tau_dag, li_dag)) return 0;

    //std::printf("\ninserting:");
    if (data.use_bare_propagator) {
      auto diagram = diagram::time_diagram_t{proposed_config.c_list, proposed_config.cdag_list, {}}; // make a free function (not member of data)
      //print_configuration(diagram);
      proposed_w.hyb  = 0; // diagram::determinant(diagram);
      proposed_u_frame = propagator_product(data.h_diag, diagram, data.tau_max);
    } else {
      auto diagram =
         diagram::time_diagram_t{proposed_config.c_list, proposed_config.cdag_list, {data.tau_split}}; // make a free function (not member of data)
      proposed_w.hyb  = 0;                                                                            //diagram::inclusion_exclusion(diagram);
      proposed_u_frame = propagator_product(data.h_diag, diagram, data.tau_max, &data.u_tau);
    }
    proposed_w.loc = frobenius_norm(proposed_u_frame);

    auto w_hyb_ratio = proposed_w.hyb / data.w.hyb;
    auto w_loc_ratio = proposed_w.loc / data.w.loc;
    auto t_ratio     = std::pow(data.tau_max / (data.size() + 1), 2);

    if (false) {
      printf("proposed_w_hyb, proposed_w_loc, last_w_hyb, last_w_loc, t_ratio:  %f %f   %f %f   %f\n", proposed_w.hyb, proposed_w.loc, data.w.hyb,
             data.w.loc, t_ratio);
      fflush(stdout);
    }
    return t_ratio * w_loc_ratio * w_hyb_ratio;
  }

  //
  scalar_t insert::accept() {
    //std::printf("yes\n");
    //print_configuration(data.get_time_diagram()  );
    data.w.hyb = proposed_w.hyb;
    data.w.loc = proposed_w.loc;

    return 1;
  }

} // namespace inchworm::moves
