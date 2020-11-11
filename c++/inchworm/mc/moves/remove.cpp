#include "./remove.hpp"

namespace inchworm::moves {

  scalar_t remove::attempt() {

    proposed_config = data.config; // we first copy last accepted config before proposing the new remove.
    proposed_w      = data.w;

    int N = data.config.size(); // size before proposition
    if (N == 0) return 0;

    int idx     = rng(N);
    int idx_dag = rng(N);

    if (not proposed_config.try_erase(idx, idx_dag)) return 0; //data is not modified in this case

    ///////

    auto diagram  = diagram::time_diagram_t{proposed_config.d_list, proposed_config.d_dag_list, {params.tau_split}};
    auto hyb_mat  = diagram::hyb_matrix_t(diagram, params.hyb_adaptor);
    proposed_sign = diagram.sign();

    if (params.use_bare_propagator)
      proposed_w.hyb = hyb_mat.det();
    else {
      if (proposed_config.size() == 0)
        proposed_w.hyb = 1;
      else
        proposed_w.hyb = diagram::inclusion_exclusion(diagram, hyb_mat);
    }

    double tol = 1e-12;
    if (std::abs(proposed_w.hyb) < tol) return 0.0;

    // FIXME: at the next big refactoring: extract into a function and all its root...
    if (params.mode == 0) {
      if (params.use_bare_propagator) {
        proposed_u_partial = impurity_product(params.ad_imp, diagram, 0, params.tau_max, nullptr);
      } else {
        // -- ip * ip
        //
        // --> proposed_green_matrix Trace( ip * d * ip * ddag )
        //
        // - <c(tau) cd(0;);>
        proposed_u_partial = impurity_product(params.ad_imp, diagram, params.tau_split, params.tau_max, &params.u_tau)
           * impurity_product(params.ad_imp, diagram, 0, params.tau_split, &params.u_tau);
      }

      proposed_w.loc = frobenius_norm(proposed_u_partial);
    } else if (params.mode == 1) {
      EXPECTS(not params.use_bare_propagator);

      auto l = impurity_product(params.ad_imp, diagram, params.tau_split, params.tau_max, &params.u_tau);
      auto r = impurity_product(params.ad_imp, diagram, 0, params.tau_split, &params.u_tau);

      for (auto &Bl : proposed_g_frame) Bl = 0;

      proposed_g_frame = make_g_frame_from_l_and_r(params.ad_imp, gf_struct, l, r);

      auto const &ops = diagram.op_list;
      int nop_r       = std::count_if(begin(ops), end(ops), [tau_split = params.tau_split](auto const &op) { return tau_split > op.tau; });
      if (nop_r % 2 == 1) {
        for (auto &bl : proposed_g_frame) bl *= -1;
      }

      proposed_w.loc = frobenius_norm(proposed_g_frame);
    }

    int n_fops       = (params.ad_imp.get_fops()).size();
    auto sign_ratio  = proposed_sign / data.sign;
    auto w_hyb_ratio = proposed_w.hyb / data.w.hyb;
    auto w_loc_ratio = proposed_w.loc / data.w.loc;
    auto t_ratio     = std::pow(N / (n_fops * params.tau_max), 2);

#ifdef INCHWORM_DEBUG_PRINTS
    std::printf("\n\n====== Try Remove ======\n");
    print_configuration(diagram);
    if (params.mode == 0)
      print(proposed_u_partial);
    else
      print(proposed_g_frame);
    hyb_mat.print();
    std::printf("\n\nhyb.det()=% 4.7f \n", hyb_mat.det());
    std::printf("\n\nsign= %d  w_hyb=% 4.7f  w_loc=% 4.7f    old_w_hyb=% 4.7f  old_w_loc=% 4.7f \n", proposed_sign, proposed_w.hyb, proposed_w.loc,
                data.w.hyb, data.w.loc);
    std::printf("\n\nsign_ratio= %d  w_hyb_ratio=% 4.7f  w_loc_ratio=% 4.7f  t_ratio=% 4.7f\n", sign_ratio, w_hyb_ratio, w_loc_ratio, t_ratio);
    if (proposed_config.size() > 0) {
      std::printf("proper_enum w.hyb         =% 4.7f \n", diagram::proper_enum(diagram, hyb_mat, 10));
      std::printf("inclusion_exclusion w.hyb =% 4.7f \n", diagram::inclusion_exclusion(diagram, hyb_mat, 10));
    }
#endif

    return sign_ratio * t_ratio * w_loc_ratio * w_hyb_ratio;
  }

  //
  scalar_t remove::accept() {
#ifdef INCHWORM_DEBUG_PRINTS
    std::printf("\n\n====== Accept Remove ======\n");
#endif
    data.w         = proposed_w;
    data.u_partial = proposed_u_partial;
    data.g_frame   = proposed_g_frame;
    data.config    = proposed_config;
    data.sign      = proposed_sign;
    return 1;
  }

} // namespace inchworm::moves
