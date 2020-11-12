#include "./insert.hpp"
#include "./../impurity_product.hpp"

namespace inchworm::moves {

  scalar_t insert::attempt() {

    prop_config  = data.config; // we first copy last accepted config before proposing the new insert
    prop_weights = data.weights;

    int n_fops = all_d_ops.size();

    auto d     = all_d_ops[rng(n_fops)];
    auto d_dag = all_d_dag_ops[rng(n_fops)];

    d.tau     = rng(params.tau_max);
    d_dag.tau = rng(params.tau_max);

    if (not prop_config.try_insert(d_dag, d)) return 0;

    ///////

    auto diagram  = diagram::time_diagram_t{prop_config, {params.tau_split}};
    auto hyb_mat  = diagram::hyb_matrix_t(diagram, params.hyb_adaptor);
    prop_sign     = diagram.sign();

    if (params.use_bare_propagator)
      prop_weights.hyb = hyb_mat.det();
    else
      prop_weights.hyb = diagram::inclusion_exclusion(diagram, hyb_mat);

    double tol = 1e-12;
    if (std::abs(prop_weights.hyb) < tol) return 0.0;

    // FIXME: at the next big refactoring: extract into a function and all its root...
    if (params.mode == 0) {
      if (params.use_bare_propagator) {
        prop_u_partial = impurity_product(params.ad_imp, diagram, 0, params.tau_max, nullptr);
      } else {
        // -- ip * ip
        //
        // --> prop_green_matrix Trace( ip * d * ip * ddag )
        //
        // - <c(tau) cd(0;);>
        prop_u_partial = impurity_product(params.ad_imp, diagram, params.tau_split, params.tau_max, &params.u_tau)
           * impurity_product(params.ad_imp, diagram, 0, params.tau_split, &params.u_tau);
      }

      prop_weights.loc = frobenius_norm(prop_u_partial);
    } else if (params.mode == 1) {
      EXPECTS(not params.use_bare_propagator);

      auto l = impurity_product(params.ad_imp, diagram, params.tau_split, params.tau_max, &params.u_tau);
      auto r = impurity_product(params.ad_imp, diagram, 0, params.tau_split, &params.u_tau);

      for (auto &Bl : prop_g_frame) Bl = 0;

      prop_g_frame = make_g_frame_from_l_and_r(params.ad_imp, gf_struct, l, r);

      auto const &ops = diagram.op_list;
      int nop_r       = std::count_if(begin(ops), end(ops), [tau_split = params.tau_split](auto const &op) { return tau_split > op.tau; });
      if (nop_r % 2 == 1) {
        for (auto &bl : prop_g_frame) bl *= -1;
      }

      prop_weights.loc = frobenius_norm(prop_g_frame);
    }

    auto sign_ratio  = prop_sign / data.sign;
    auto w_hyb_ratio = prop_weights.hyb / data.weights.hyb;
    auto w_loc_ratio = prop_weights.loc / data.weights.loc;

    int N        = data.config.size(); // size before proposition
    auto t_ratio = std::pow(params.tau_max * n_fops / (N + 1), 2);

#ifdef INCHWORM_DEBUG_PRINTS
    std::printf("\n\n====== Try Insert ======\n");
    print_configuration(diagram);
    if (params.mode == 0)
      print(prop_u_partial);
    else
      print(prop_g_frame);
    hyb_mat.print();
    std::printf("\n\nhyb.det()=% 4.7f \n", hyb_mat.det());
    std::printf("\n\nsign= %d  w_hyb=% 4.7f  w_loc=% 4.7f    old_w_hyb=% 4.7f  old_w_loc=% 4.7f \n", prop_sign, prop_weights.hyb, prop_weights.loc,
                data.w.hyb, data.w.loc);
    std::printf("\n\nsign_ratio= %d  w_hyb_ratio=% 4.7f  w_loc_ratio=% 4.7f  t_ratio=% 4.7f\n", sign_ratio, w_hyb_ratio, w_loc_ratio, t_ratio);
    std::printf("proper_enum w.hyb         =% 4.7f \n", diagram::proper_enum(diagram, hyb_mat, 10));
    std::printf("inclusion_exclusion w.hyb =% 4.7f \n", diagram::inclusion_exclusion(diagram, hyb_mat, 10));
#endif

    return sign_ratio * t_ratio * w_loc_ratio * w_hyb_ratio;
  }

  //
  scalar_t insert::accept() {
#ifdef INCHWORM_DEBUG_PRINTS
    std::printf("\n\n====== Accept Insert ======\n");
#endif
    data.weights   = prop_weights;
    data.u_partial = prop_u_partial;
    data.g_frame   = prop_g_frame;
    data.config    = prop_config;
    data.sign      = prop_sign;
    return 1;
  }

} // namespace inchworm::moves
