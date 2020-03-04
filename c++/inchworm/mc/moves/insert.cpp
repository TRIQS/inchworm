#include "./insert.hpp"

namespace inchworm::moves {

  scalar_t insert::attempt() {

    proposed_config = data.config; // we first copy last accepted config before proposing the new insert
    proposed_w      = data.w;

    int N      = data.config.size(); // size before proposition
    int n_fops = (params.h_diag.get_fops()).size();
    int li     = rng(n_fops);
    int li_dag = rng(n_fops);

    double tau     = rng(params.tau_max);
    double tau_dag = rng(params.tau_max);
    if (not proposed_config.try_insert(tau, li, tau_dag, li_dag)) return 0;

    auto diagram  = diagram::time_diagram_t{proposed_config.c_list, proposed_config.cdag_list, {params.tau_split}};
    auto hyb_mat  = diagram::hyb_matrix_t(diagram, params.hyb_adaptor);
    proposed_sign = diagram.sign();

    if (params.use_bare_propagator)
      proposed_w.hyb = hyb_mat.det();
    else
      proposed_w.hyb = diagram::inclusion_exclusion(diagram, hyb_mat);

    double tol = 1e-12;
    if (std::abs(proposed_w.hyb) < tol) return 0.0;
    //if (proposed_config.size() >5) return 0.0;

    if (params.use_bare_propagator)
      proposed_u_frame = propagator_product(params.h_diag, diagram, params.tau_max);
    else
      proposed_u_frame = propagator_product(params.h_diag, diagram, params.tau_max, &params.u_tau);

    proposed_w.loc = frobenius_norm(proposed_u_frame);

    proposed_w.loc   = frobenius_norm(proposed_u_frame);
    auto sign_ratio  = proposed_sign / data.sign;
    auto w_hyb_ratio = proposed_w.hyb / data.w.hyb;
    auto w_loc_ratio = proposed_w.loc / data.w.loc;
    auto t_ratio     = std::pow(params.tau_max * n_fops / (N + 1), 2);

    /*if (proposed_config.size() == 1) {
      std::printf("\n ");
      for (auto &B : proposed_u_frame) { std::cout << B; }
      hyb_mat.print();
      std::printf("\n\nhyb.det()=% 4.7f \n", hyb_mat.det()),
         std::printf("\n\nsign= %d  w_hyb=% 4.7f  w_loc=% 4.7f    old_w_hyb=% 4.7f  old_w_loc=% 4.7f \n", proposed_sign, proposed_w.hyb,
                     proposed_w.loc, data.w.hyb, data.w.loc);
      std::printf("\n\nsign_ratio= %d  w_hyb_ratio=% 4.7f  w_loc_ratio=% 4.7f  t_ratio=% 4.7f\n", sign_ratio, w_hyb_ratio, w_loc_ratio, t_ratio);
    }*/
    return sign_ratio * t_ratio * w_loc_ratio * w_hyb_ratio;
  }

  //
  scalar_t insert::accept() {
    //std::printf("\n\nsize=%d\n", proposed_config.size());
    //for (auto const &B : proposed_u_frame) std::cout << B;
    if (proposed_config.size() == 100) { //params.verbosity == 10) {
      auto diagram = diagram::time_diagram_t(proposed_config.c_list, proposed_config.cdag_list, {});
      print_configuration(diagram);
      auto hyb_mat = diagram::hyb_matrix_t(diagram, params.hyb_adaptor);
      hyb_mat.print();
      //std::printf("\n\nsign= %d  t_ratio=% 4.7f  w_hyb=% 4.7f  w_loc=% 4.7f\n", sign, params.tau_max / (proposed_config.size() + 1), proposed_w.hyb, proposed_w.loc);
      for (auto const &B : proposed_u_frame) std::cout << B;
      std::printf("\n\nsign= % d   w_hyb=% 4.7f  w_loc=% 4.7f\n\n\n", proposed_sign, proposed_w.hyb, proposed_w.loc);
      std::printf("\n\nhyb.det()=% 4.7f \n", hyb_mat.det());
    }
    data.w       = proposed_w;
    data.u_frame = proposed_u_frame;
    data.config  = proposed_config;
    data.sign    = proposed_sign;
    //std::swap(proposed_w, data.w);
    //std::swap(proposed_u_frame, data.u_frame);
    //std::swap(proposed_config, data.config);
    return 1;
  }

} // namespace inchworm::moves
