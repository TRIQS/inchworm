#include "./remove.hpp"
#include "../diagram/print.hpp"

namespace inchworm::moves {

  scalar_t remove::attempt() {

    // temporary values:
    proposed_config = data.config; // we first copy last accepted config before proposing the new remove.
    proposed_w      = data.w;

    int N     = data.size();
    int i     = rng(N);
    int i_dag = rng(N);
    if (not proposed_config.try_erase(i, i_dag)) return 0; //data is not modified in this case

    //std::printf("success? %d \n", success);
    //}

    //std::printf("\nremoving:");
    if (data.use_bare_propagator) {
      auto diagram = diagram::time_diagram_t{proposed_config.c_list, proposed_config.cdag_list, {}}; // make a free function (not member of data)
      //print_configuration(diagram);
      proposed_w.hyb   = 0; // diagram::determinant(diagram);
      proposed_u_frame = propagator_product(data.h_diag, diagram, data.tau_max);
    } else {
      auto diagram =
         diagram::time_diagram_t{proposed_config.c_list, proposed_config.cdag_list, {data.tau_split}}; // make a free function (not member of data)
      proposed_w.hyb   = 0;                                                                            //diagram::inclusion_exclusion(diagram);
      proposed_u_frame = propagator_product(data.h_diag, diagram, data.tau_max, &data.u_tau);
    }
    proposed_w.loc = frobenius_norm(proposed_u_frame);

    auto w_hyb_ratio = proposed_w.hyb / data.w.hyb;
    auto w_loc_ratio = proposed_w.loc / data.w.loc;
    auto t_ratio     = std::pow(N / data.tau_max, 2);

    //printf("helloaaaaaaaaaaiii02:  %f  \n", t_ratio * w_loc_ratio * w_hyb_ratio);
    if (false) {
      printf("proposed_w_hyb, proposed_w_loc, last_w_hyb, last_w_loc, t_ratio:  %f %f   %f %f   %f\n", proposed_w.hyb, proposed_w.loc, data.w.hyb,
             data.w.loc, t_ratio);
      fflush(stdout);
    }
    return t_ratio * w_loc_ratio * w_hyb_ratio;
  }

  //
  scalar_t remove::accept() {
    //std::printf("yes\n");
    //print_configuration(data.get_time_diagram()  );
    data.w.hyb = proposed_w.hyb;
    data.w.loc = proposed_w.loc;

    return 1;
  }

} // namespace inchworm::moves
