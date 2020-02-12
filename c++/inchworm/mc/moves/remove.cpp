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

    int sign;
    if (params.use_bare_propagator) {
      auto diagram = diagram::time_diagram_t(proposed_config.c_list, proposed_config.cdag_list, {}); // make a free function (not member of data)
      //print_configuration(diagram);
      auto hyb_mat     = diagram::hyb_matrix_t(diagram, params.hyb_adaptor);
      //std::printf("\n\nsize=%d,   det=% 4.5e\n", proposed_config.size(), hyb_mat.det());
      
      proposed_w.hyb   = hyb_mat.det();
      proposed_u_frame = propagator_product(params.h_diag, diagram, params.tau_max);
      sign             = diagram.sign();
    } else {
      auto diagram =
         diagram::time_diagram_t{proposed_config.c_list, proposed_config.cdag_list, {params.tau_split}}; // make a free function (not member of data)
      auto hyb_mat     = diagram::hyb_matrix_t(diagram, params.hyb_adaptor);
      proposed_w.hyb   = diagram::inclusion_exclusion(diagram, hyb_mat);
      proposed_u_frame = propagator_product(params.h_diag, diagram, params.tau_max, &params.u_tau);
      sign             = diagram.sign();
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
    return sign * t_ratio * w_loc_ratio * w_hyb_ratio;
  }

  //
  scalar_t remove::accept() {
    //std::printf("yes\n");
    //std::printf("\n\nsize=%d\n", proposed_config.size());
    //for (auto const &B : proposed_u_frame) std::cout << B;
    data.w       = proposed_w;
    data.u_frame = proposed_u_frame;
    data.config  = proposed_config;
    return 1;
  }

} // namespace inchworm::moves
