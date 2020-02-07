#include "./insert.hpp"

namespace inchworm::moves {

  scalar_t insert::attempt() {

    proposed_config = data.config; // we first copy last accepted config before proposing the new insert
    proposed_w      = data.w;

    int N      = proposed_config.size(); // size before proposition
    int n_fops = (params.h_diag.get_fops()).size();
    int li     = rng(n_fops);
    int li_dag = rng(n_fops);

    double tau     = rng(params.tau_max);
    double tau_dag = rng(params.tau_max);
    if (not proposed_config.try_insert(tau, li, tau_dag, li_dag)) return 0;

    //std::printf("\ninserting:");
    if (params.use_bare_propagator) {
      auto diagram = diagram::time_diagram_t{proposed_config.c_list, proposed_config.cdag_list, {}}; // make a free function (not member of data)
      auto hyb_mat = diagram::hyb_matrix_t(diagram);
      //print_configuration(diagram);
      proposed_w.hyb   = hyb_mat.det();
      proposed_u_frame = propagator_product(params.h_diag, diagram, params.tau_max);
    } else {
      auto diagram =
         diagram::time_diagram_t{proposed_config.c_list, proposed_config.cdag_list, {params.tau_split}}; // make a free function (not member of data)
      auto hyb_mat     = diagram::hyb_matrix_t(diagram);
      proposed_w.hyb   = diagram::inclusion_exclusion(diagram, hyb_mat);
      proposed_u_frame = propagator_product(params.h_diag, diagram, params.tau_max, &params.u_tau);
    }
    proposed_w.loc = frobenius_norm(proposed_u_frame);

    auto w_hyb_ratio = proposed_w.hyb / data.w.hyb;
    auto w_loc_ratio = proposed_w.loc / data.w.loc;
    auto t_ratio     = std::pow(params.tau_max / (N + 1), 2);

    if (false) {
      std::printf("proposed_w_hyb, proposed_w_loc, last_w_hyb, last_w_loc, t_ratio:  %f %f   %f %f   %f\n", proposed_w.hyb, proposed_w.loc,
                  data.w.hyb, data.w.loc, t_ratio);
      fflush(stdout);
    }
    return t_ratio * w_loc_ratio * w_hyb_ratio;
  }

  //
  scalar_t insert::accept() {
    //std::printf("yes\n");
    //print_configuration(data.get_time_diagram()  );
    data.w       = proposed_w;
    data.u_frame = proposed_u_frame;
    data.config  = proposed_config;
    //std::swap(proposed_w, data.w);
    //std::swap(proposed_u_frame, data.u_frame);
    //std::swap(proposed_config, data.config);
    return 1;
  }

} // namespace inchworm::moves
