#include "./double_insert.hpp"

namespace inchworm::moves {

  scalar_t double_insert::attempt() {

    proposed_config = data.config; // we first copy last accepted config before proposing the new double_insert
    proposed_w      = data.w;

    int N      = data.config.size(); // size before proposition
    int n_fops = (params.ad_imp.get_fops()).size();

    int li1     = rng(n_fops);
    int li1_dag = rng(n_fops);
    int li2     = rng(n_fops);
    int li2_dag = rng(n_fops);

    double tau1     = rng(params.tau_max);
    double tau1_dag = rng(params.tau_max);
    double tau2     = rng(params.tau_max);
    double tau2_dag = rng(params.tau_max);

    if (not proposed_config.try_double_insert(tau1, li1, tau1_dag, li1_dag, tau2, li2, tau2_dag, li2_dag)) return 0;

    ///////

    auto diagram  = diagram::time_diagram_t{proposed_config.d_list, proposed_config.d_dag_list, {params.tau_split}};
    auto hyb_mat  = diagram::hyb_matrix_t(diagram, params.hyb_adaptor);
    proposed_sign = diagram.sign();

    if (params.use_bare_propagator)
      proposed_w.hyb = hyb_mat.det();
    else
      //proposed_w.hyb = diagram::proper_enum(diagram, hyb_mat);
      proposed_w.hyb = diagram::inclusion_exclusion(diagram, hyb_mat);

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
      //proposed_w.loc = trace(proposed_u_partial);

    } else if (params.mode == 1) {
      EXPECTS(not params.use_bare_propagator);

      auto l = impurity_product(params.ad_imp, diagram, params.tau_split, params.tau_max, &params.u_tau);
      auto r = impurity_product(params.ad_imp, diagram, 0, params.tau_split, &params.u_tau);

      for (auto &Bl : proposed_g_frame) Bl = 0;

      proposed_g_frame = make_g_frame_from_l_and_r(params.ad_imp, params.map_lin_idx_to_block_inner, gf_struct, l, r);

      auto const &ops = diagram.op_list;
      int nop_r       = std::count_if(begin(ops), end(ops), [tau_split = params.tau_split](auto const &op) { return tau_split > op.tau; });
      if (nop_r % 2 == 1) {
        for (auto &bl : proposed_g_frame) bl *= -1;
      }

      proposed_w.loc = frobenius_norm(proposed_g_frame);
    }

    auto sign_ratio  = proposed_sign / data.sign;
    auto w_hyb_ratio = proposed_w.hyb / data.w.hyb;
    auto w_loc_ratio = proposed_w.loc / data.w.loc;
    auto t_ratio     = std::pow(params.tau_max * n_fops / (N + 2), 4);

#ifdef INCHWORM_DEBUG_PRINTS
    std::printf("\n\n====== Try Insert2 ======\n");
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
    std::printf("proper_enum w.hyb         =% 4.7f \n", diagram::proper_enum(diagram, hyb_mat, 10));
    std::printf("inclusion_exclusion w.hyb =% 4.7f \n", diagram::inclusion_exclusion(diagram, hyb_mat, 10));
#endif

    return sign_ratio * t_ratio * w_loc_ratio * w_hyb_ratio;
  }

  //
  scalar_t double_insert::accept() {
#ifdef INCHWORM_DEBUG_PRINTS
    std::printf("\n\n====== Accept Insert2 ======\n");
#endif
    data.w         = proposed_w;
    data.u_partial = proposed_u_partial;
    data.g_frame   = proposed_g_frame;
    data.config    = proposed_config;
    data.sign      = proposed_sign;
    return 1;
  }

} // namespace inchworm::moves
