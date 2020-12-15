#include "./base_move.hpp"
#include "./../diagram/hyb_matrix.hpp"
#include "./../diagram/inclusion_exclusion.hpp"
#include "./../u_frame.hpp"
#include "./../atom_diag.hpp"
#include "./../impurity_product.hpp"

namespace inchworm::moves {

  base_move::base_move(qmc_config_data_t &data, gf_struct_t const &gf_struct, qmc_params_t const &qmc_params, triqs::mc_tools::random_generator &rng)
     : data(data), prop_data(data), params(qmc_params), rng(rng), gf_struct(gf_struct) {

    for (auto const &op : qmc_params.ad_imp.get_fops()) {
      auto bl_name = std::get<std::string>(op.index[0]);
      auto idx     = std::get<long>(op.index[1]);

      // Determine the number of the bl_name in gf_struct
      auto it = std::find_if(gf_struct.cbegin(), gf_struct.cend(), [&](auto &&x) { return x.first == bl_name; });
      long bl = std::distance(gf_struct.cbegin(), it);

      all_d_ops.push_back({0.0, false, op.linear_index, bl, idx});
      all_d_dag_ops.push_back({0.0, true, op.linear_index, bl, idx});
    }
  }

  scalar_t base_move::attempt() {

    prop_data = data;

    // ------ Generate the new configuration -------

    auto t_ratio = try_config_update(prop_data.config);
    if (t_ratio == 0.0) return 0.0; // Check if move has failed

    // ------ Calculate the hybridization weight -------

    auto diagram   = diagram::time_diagram_t{prop_data.config, {params.tau_split}};
    auto hyb_mat   = diagram::hyb_matrix_t(diagram, params.hyb_adaptor);
    prop_data.sign = diagram.sign();

    if (params.use_bare_propagator)
      prop_data.weights.hyb = hyb_mat.det();
    else
      prop_data.weights.hyb = diagram::inclusion_exclusion(diagram, hyb_mat);

    double tol = 1e-12;
    if (std::abs(prop_data.weights.hyb) < tol) return 0.0;

    // ------ Calculate the impurity frame weight -------

    if (params.mode == MODE::PROPAGATOR) {
      if (params.use_bare_propagator) {
        prop_data.frame = make_frame(impurity_product(params.ad_imp, diagram, params.tau_max, 0));
      } else { // FIXME incorporate treatment of tau_split into impurity product
        prop_data.frame = make_frame(impurity_product(params.ad_imp, diagram, params.tau_max, params.tau_split, &params.u_tau)
                                     * impurity_product(params.ad_imp, diagram, params.tau_split, 0, &params.u_tau));
      }

    } else { // MODE::GREENFUNCTION
      EXPECTS(not params.use_bare_propagator);

      prop_data.frame = make_zero_frame(gf_struct);

      // Calculate -Tr[imp_prod(beta, tau) * c(tau) * imp_prod(tau, 0) * cdag(0)]
      // for all combinations of fundamental operator flavors
      auto l          = impurity_product(params.ad_imp, diagram, params.tau_max, params.tau_split, &params.u_tau);
      auto r          = impurity_product(params.ad_imp, diagram, params.tau_split, 0, &params.u_tau);
      prop_data.frame = make_g_frame_from_l_and_r(params.ad_imp, gf_struct, l, r);

      // Account for the sign due to the additional operator insertions
      auto const &ops = diagram.op_list;
      int nop_r       = std::count_if(begin(ops), end(ops), [tau_split = params.tau_split](auto const &op) { return tau_split > op.tau; });
      if (nop_r % 2 == 1) { prop_data.frame *= -1; }
    }

    prop_data.weights.imp = frobenius_norm(prop_data.frame);

    // Reweight perturbation orders below the reweighting_cutoff to guarantee
    // that the zeroth order is sampled properly for normalization purposes
    // FIXME Improve reweighting using proper perturbation order histogram
    if (diagram.perturbation_order() < reweighting_cutoff)
      return prop_data.weights.imp *= std::pow(reweighting_cutoff - diagram.perturbation_order(), 2);

    // ------ Calculate overall weight ratio -------

    auto sign_ratio  = prop_data.sign / data.sign;
    auto w_hyb_ratio = prop_data.weights.hyb / data.weights.hyb;
    auto w_imp_ratio = prop_data.weights.imp / data.weights.imp;

    auto ratio = sign_ratio * t_ratio * w_imp_ratio * w_hyb_ratio;

    // ------ Debugging Information -------

#ifdef INCHWORM_DEBUG_PRINTS
    std::printf("\n\n====== Try %s ======\n", name());
    print_configuration(diagram);
    if (params.mode == MODE::PROPAGATOR)
      print(prop_data.u_partial);
    else // MODE::GREENFUNCTION
      print(prop_data.g_frame);
    hyb_mat.print();
    std::printf("\n\nhyb.det()=% 4.7f \n", hyb_mat.det());
    std::printf("\n\nsign= %d  w_hyb=% 4.7f  w_imp=% 4.7f    old_w_hyb=% 4.7f  old_w_imp=% 4.7f \n", prop_data.sign, prop_data.weights.hyb,
                prop_data.weights.imp, data.w.hyb, data.w.imp);
    std::printf("\n\nsign_ratio= %d  w_hyb_ratio=% 4.7f  w_imp_ratio=% 4.7f  t_ratio=% 4.7f\n", sign_ratio, w_hyb_ratio, w_imp_ratio, t_ratio);
    std::printf("proper_enum w.hyb         =% 4.7f \n", diagram::proper_enum(diagram, hyb_mat, true /*verbose*/));
    std::printf("inclusion_exclusion w.hyb =% 4.7f \n", diagram::inclusion_exclusion(diagram, hyb_mat, true /*verbose*/));
#endif

    return ratio;
  }

  scalar_t base_move::accept() {
#ifdef INCHWORM_DEBUG_PRINTS
    std::printf("\n\n====== Accept %s ======\n", name());
#endif
    data = prop_data;
    return 1.0;
  }

} // namespace inchworm::moves
