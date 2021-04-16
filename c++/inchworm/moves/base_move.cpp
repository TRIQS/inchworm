#include "./base_move.hpp"

#include "./../diagram/hyb_matrix.hpp"
#include "./../diagram/proper_enum.hpp"
#include "./../diagram/inclusion_exclusion.hpp"
#include "./../diagram/print.hpp"
#include "./../u_frame.hpp"
#include "./../atom_diag.hpp"
#include "./../impurity_product.hpp"
#include "./../util.hpp"

namespace inchworm::moves {

  scalar_t base_move::attempt() {

    prop_config = config;

    // ------ Generate the new configuration and diagram -------

    auto t_ratio = try_config_update(prop_config);
    if (t_ratio == 0.0) return 0.0; // Check if move has failed
    auto prop_pert_order = prop_config.size();
    if (prop_pert_order > solver.last_solve_params->max_order.value_or(prop_pert_order)) return 0.0;
    auto diagram = diagram::time_diagram_t{prop_config, {params.tau_split}};

    // We need to have at least one split-point between operators for a finite hybridization weight
    if (not params.use_bare_propagator and diagram.size() > 0 and diagram.is_trivial) return 0.0;

    // Quick-check for vanishing impurity trace
    if (params.mode == MODE::PROPAGATOR and has_zero_trace(solver.ad_imp, diagram)) return 0.0;

    // ------ Calculate the hybridization weight -------

    auto hyb_mat     = diagram::hyb_matrix_t(diagram, solver.Delta_tau);
    prop_config.sign = diagram.sign();

    if (params.use_bare_propagator)
      prop_config.hyb_weight = hyb_mat.det();
    else
#ifdef PROPER_ENUMERATION
      prop_config.hyb_weight = diagram::proper_enum(diagram, hyb_mat);
#else
      prop_config.hyb_weight = diagram::inclusion_exclusion(diagram, hyb_mat);
#endif

    double tol = 1e-12;
    if (std::abs(prop_config.hyb_weight) < tol) return 0.0;

    // ------ Calculate the impurity frame weight -------

    prop_frame = frame;
    if (params.mode == MODE::PROPAGATOR) {
      if (params.use_bare_propagator) {
        prop_frame = make_frame(impurity_product(solver.ad_imp, diagram, params.tau_max, 0));
      } else { // FIXME incorporate treatment of tau_split into impurity product
        prop_frame = make_frame(impurity_product(solver.ad_imp, diagram, params.tau_max, params.tau_split, &solver.u_interpolator)
                                * impurity_product(solver.ad_imp, diagram, params.tau_split, 0, &solver.u_interpolator));
      }

    } else { // MODE::GREENFUNCTION
      EXPECTS(not params.use_bare_propagator);

      prop_frame = make_zero_frame(gf_struct);

      // Calculate -Tr[imp_prod(beta, tau) * c(tau) * imp_prod(tau, 0) * cdag(0)]
      // for all combinations of fundamental operator flavors
      auto l          = impurity_product(solver.ad_imp, diagram, params.tau_max, params.tau_split, &solver.u_interpolator);
      auto r          = impurity_product(solver.ad_imp, diagram, params.tau_split, 0, &solver.u_interpolator);
      prop_frame      = make_g_frame_from_l_and_r(solver.ad_imp, gf_struct, l, r);

      // Account for the sign due to the additional operator insertions
      auto const &ops = diagram.op_list;
      int nop_r       = std::count_if(begin(ops), end(ops), [tau_split = params.tau_split](auto const &op) { return tau_split > op.tau; });
      if (nop_r % 2 == 1) { prop_frame *= -1; }
    }

    prop_config.imp_weight = frobenius_norm(prop_frame);

    // Reweight perturbation orders below the reweighting_cutoff to guarantee
    // that the zeroth order is sampled properly for normalization purposes
    if (prop_pert_order < reweighting_cutoff) prop_config.imp_weight *= reweighting_coeffs[prop_pert_order];

    // ------ Calculate overall weight ratio -------

    auto sign_ratio  = prop_config.sign / config.sign;
    auto w_hyb_ratio = prop_config.hyb_weight / config.hyb_weight;
    auto w_imp_ratio = prop_config.imp_weight / config.imp_weight;

    auto ratio = sign_ratio * t_ratio * w_imp_ratio * w_hyb_ratio;

    // ------ Debugging Information -------

#ifdef INCHWORM_DEBUG_PRINTS
    //if(rng.preview() >= std::min(1.0, std::abs(ratio))) return ratio;
    std::printf("\n\n====== Try %s ======\n", name().c_str());
    print_configuration(diagram);
    print(prop_frame);
    hyb_mat.print();
    std::printf("\n\nhyb.det()=% 4.7e \n", hyb_mat.det());
    std::printf("\n\nsign= %d  w_hyb=% 4.7e  w_imp=% 4.7e    old_w_hyb=% 4.7e  old_w_imp=% 4.7e \n", prop_config.sign, prop_config.hyb_weight,
                prop_config.imp_weight, config.hyb_weight, config.imp_weight);
    std::printf("\n\nsign_ratio= %d  w_hyb_ratio=% 4.7e  w_imp_ratio=% 4.7e  t_ratio=% 4.7e\n", sign_ratio, w_hyb_ratio, w_imp_ratio, t_ratio);
    std::printf("\n\nratio=% 4.7e\n", ratio);
    if (not params.use_bare_propagator) {
      std::printf("proper_enum w.hyb         =% 4.7f \n", diagram::proper_enum(diagram, hyb_mat));
      std::printf("inclusion_exclusion w.hyb =% 4.7f \n", diagram::inclusion_exclusion(diagram, hyb_mat));
    }
    getchar();
#endif

    return ratio;
  }

  scalar_t base_move::accept() {
#ifdef INCHWORM_DEBUG_PRINTS
    std::printf("\n\n====== Accept %s ======\n", name().c_str());
#endif
    config = prop_config;
    frame  = prop_frame;
    return 1.0;
  }

} // namespace inchworm::moves
