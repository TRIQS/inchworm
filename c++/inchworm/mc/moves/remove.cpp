#include "./remove.hpp"

namespace inchworm::moves {

  scalar_t remove::attempt() {

    proposed_config = data.config; // we first copy last accepted config before proposing the new remove.
    proposed_w      = data.w;

    int N     = proposed_config.size(); // size before proposition
    int i     = rng(N);
    int i_dag = rng(N);
    if (not proposed_config.try_erase(i, i_dag)) return 0; //data is not modified in this case

    //std::printf("\nremoving:");
    if (params.use_bare_propagator) {
      auto diagram = diagram::time_diagram_t{proposed_config.c_list, proposed_config.cdag_list, {}}; // make a free function (not member of data)
      //print_configuration(diagram);
      proposed_w.hyb   = 0; // diagram::determinant(diagram);
      proposed_u_frame = propagator_product(params.h_diag, diagram, params.tau_max);
    } else {
      auto diagram =
         diagram::time_diagram_t{proposed_config.c_list, proposed_config.cdag_list, {params.tau_split}}; // make a free function (not member of data)
      proposed_w.hyb   = 0;                                                                              //diagram::inclusion_exclusion(diagram);
      proposed_u_frame = propagator_product(params.h_diag, diagram, params.tau_max, &params.u_tau);
    }
    proposed_w.loc = frobenius_norm(proposed_u_frame);

    auto w_hyb_ratio = proposed_w.hyb / data.w.hyb;
    auto w_loc_ratio = proposed_w.loc / data.w.loc;
    auto t_ratio     = std::pow(N / params.tau_max, 2);

    if (false) {
      std::printf("proposed_w_hyb, proposed_w_loc, last_w_hyb, last_w_loc, t_ratio:  %f %f   %f %f   %f\n", proposed_w.hyb, proposed_w.loc,
                  data.w.hyb, data.w.loc, t_ratio);
      fflush(stdout);
    }
    return t_ratio * w_loc_ratio * w_hyb_ratio;
  }

  //
  scalar_t remove::accept() {
    //std::printf("yes\n");
    //print_configuration(data.get_time_diagram()  );
    std::swap(data.w, proposed_w);
    std::swap(proposed_u_frame, data.u_frame);
    return 1;
  }

} // namespace inchworm::moves
