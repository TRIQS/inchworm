#include "./remove.hpp"

namespace inchworm::moves {

  scalar_t remove::attempt() {

    proposed_config = data.config; // we first copy last accepted config before proposing the new remove.
    proposed_w      = data.w;

    int N = data.config.size(); // size before proposition
    if (N == 0) return 0;
    int i     = rng(N);
    int i_dag = rng(N);
    //std::printf("i=%d i_dag=%d N=%d  ", i, i_dag, N);
    if (not proposed_config.try_erase(i, i_dag)) return 0; //data is not modified in this case

    //std::printf("\nremoving:");
    auto diagram  = diagram::time_diagram_t{proposed_config.c_list, proposed_config.cdag_list, {params.tau_split}};
    auto hyb_mat  = diagram::hyb_matrix_t(diagram, params.hyb_adaptor);
    proposed_sign = diagram.sign();

    if (params.use_bare_propagator)
      proposed_w.hyb = hyb_mat.det();
    else
      proposed_w.hyb = diagram::inclusion_exclusion(diagram, hyb_mat);

    double tol = 1e-12;
    if (std::abs(proposed_w.hyb) < tol) return 0.0;

    if (params.use_bare_propagator)
      proposed_u_frame = propagator_product(params.h_diag, diagram, params.tau_max);
    else
      proposed_u_frame = propagator_product(params.h_diag, diagram, params.tau_max, &params.u_tau);

    proposed_w.loc = frobenius_norm(proposed_u_frame);

    int n_fops       = (params.h_diag.get_fops()).size();
    auto sign_ratio  = proposed_sign / data.sign;
    auto w_hyb_ratio = proposed_w.hyb / data.w.hyb;
    auto w_loc_ratio = proposed_w.loc / data.w.loc;
    auto t_ratio     = std::pow(N / (n_fops * params.tau_max), 2);

    //std::printf("\n\nsign_ratio= %d  w_hyb_ratio=% 4.7f  w_loc_ratio=% 4.7f  t_ratio=% 4.7f\n", sign_ratio, w_hyb_ratio, w_loc_ratio, t_ratio);
    return sign_ratio * t_ratio * w_loc_ratio * w_hyb_ratio;
  }

  //
  scalar_t remove::accept() {
    //std::printf("yes\n");
    //std::printf("\n\nsize=%d\n", proposed_config.size());
    //for (auto const &B : proposed_u_frame) std::cout << B;
    if (proposed_config.size()==100) { //params.verbosity == 10) {
      auto diagram = diagram::time_diagram_t(proposed_config.c_list, proposed_config.cdag_list, {});
      print_configuration(diagram);
      auto hyb_mat = diagram::hyb_matrix_t(diagram, params.hyb_adaptor);
      hyb_mat.print();
      //std::printf("\n\nsign= %d  t_ratio=% 4.7f  w_hyb=% 4.7f  w_loc=% 4.7f\n", sign, params.tau_max / (proposed_config.size() + 1), proposed_w.hyb, proposed_w.loc);
      for (auto const &B : proposed_u_frame) std::cout << B;
      std::printf("\n\nsign= % d   w_hyb=% 4.7f  w_loc=% 4.7f\n\n\n", proposed_sign, proposed_w.hyb, proposed_w.loc);
    }
    data.w       = proposed_w;
    data.u_frame = proposed_u_frame;
    data.config  = proposed_config;
    data.sign    = proposed_sign;
    return 1;
  }

} // namespace inchworm::moves
