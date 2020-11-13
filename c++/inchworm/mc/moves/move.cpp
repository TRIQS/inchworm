#include "./move.hpp"
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

    auto t_ratio = try_move(prop_data.config);
    if (t_ratio == 0.0) return 0.0;

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

    if (params.mode == 0) { // --- Propagator Mode
      if (params.use_bare_propagator) {
        prop_data.u_partial = impurity_product(params.ad_imp, diagram, 0, params.tau_max, nullptr);
      } else {
        // -- ip * ip
        //
        // --> prop_data.green_matrix Trace( ip * d * ip * ddag )
        //
        // - <c(tau) cd(0;);>
        prop_data.u_partial = impurity_product(params.ad_imp, diagram, params.tau_split, params.tau_max, &params.u_tau)
           * impurity_product(params.ad_imp, diagram, 0, params.tau_split, &params.u_tau);
      }

      prop_data.weights.loc = frobenius_norm(prop_data.u_partial);

    } else if (params.mode == 1) { // --- Green Function Mode
      EXPECTS(not params.use_bare_propagator);

      prop_data.g_frame = make_frame(gf_struct);

      auto l = impurity_product(params.ad_imp, diagram, params.tau_split, params.tau_max, &params.u_tau);
      auto r = impurity_product(params.ad_imp, diagram, 0, params.tau_split, &params.u_tau);

      prop_data.g_frame = make_g_frame_from_l_and_r(params.ad_imp, gf_struct, l, r);

      // Account for the sign due to additional operator insertions
      auto const &ops = diagram.op_list;
      int nop_r       = std::count_if(begin(ops), end(ops), [tau_split = params.tau_split](auto const &op) { return tau_split > op.tau; });
      if (nop_r % 2 == 1) {
        for (auto &bl : prop_data.g_frame) bl *= -1;
      }

      prop_data.weights.loc = frobenius_norm(prop_data.g_frame);
    }

    // ------ Calculate overall weight ratio -------

    auto sign_ratio  = prop_data.sign / data.sign;
    auto w_hyb_ratio = prop_data.weights.hyb / data.weights.hyb;
    auto w_loc_ratio = prop_data.weights.loc / data.weights.loc;

    auto ratio = sign_ratio * t_ratio * w_loc_ratio * w_hyb_ratio;

    // ------ Debugging Information -------

#ifdef INCHWORM_DEBUG_PRINTS
    std::printf("\n\n====== Try %s ======\n", name());
    print_configuration(diagram);
    if (params.mode == 0)
      print(prop_data.u_partial);
    else
      print(prop_data.g_frame);
    hyb_mat.print();
    std::printf("\n\nhyb.det()=% 4.7f \n", hyb_mat.det());
    std::printf("\n\nsign= %d  w_hyb=% 4.7f  w_loc=% 4.7f    old_w_hyb=% 4.7f  old_w_loc=% 4.7f \n", prop_data.sign, prop_data.weights.hyb,
                prop_data.weights.loc, data.w.hyb, data.w.loc);
    std::printf("\n\nsign_ratio= %d  w_hyb_ratio=% 4.7f  w_loc_ratio=% 4.7f  t_ratio=% 4.7f\n", sign_ratio, w_hyb_ratio, w_loc_ratio, t_ratio);
    std::printf("proper_enum w.hyb         =% 4.7f \n", diagram::proper_enum(diagram, hyb_mat, 10));
    std::printf("inclusion_exclusion w.hyb =% 4.7f \n", diagram::inclusion_exclusion(diagram, hyb_mat, 10));
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

  // --------

  scalar_t insert::try_move(config_t &config) {
    int n_fops = all_d_ops.size();
    auto d     = all_d_ops[rng(n_fops)];
    auto d_dag = all_d_dag_ops[rng(n_fops)];

    d.tau     = rng(params.tau_max);
    d_dag.tau = rng(params.tau_max);

    if (not config.try_insert(d_dag, d)) return 0.0;

    int N = config.size();
    return std::pow(params.tau_max * n_fops / N, 2);
  }

  scalar_t double_insert::try_move(config_t &config) {
    int n_fops  = all_d_ops.size();
    auto d1     = all_d_ops[rng(n_fops)];
    auto d1_dag = all_d_dag_ops[rng(n_fops)];
    auto d2     = all_d_ops[rng(n_fops)];
    auto d2_dag = all_d_dag_ops[rng(n_fops)];

    d1.tau     = rng(params.tau_max);
    d1_dag.tau = rng(params.tau_max);
    d2.tau     = rng(params.tau_max);
    d2_dag.tau = rng(params.tau_max);

    if (not config.try_double_insert(d1_dag, d1, d2_dag, d2)) return 0;

    int N = config.size();
    return std::pow(params.tau_max * n_fops / N, 4);
  }

  scalar_t remove::try_move(config_t &config) {
    int N = config.size();
    if (N == 0) return 0;

    int idx     = rng(N);
    int idx_dag = rng(N);

    if (not config.try_erase(idx, idx_dag)) return 0; //data is not modified in this case

    int n_fops = all_d_ops.size();
    return std::pow(N / (n_fops * params.tau_max), 2);
  }

  scalar_t double_remove::try_move(config_t &config) {
    int N = config.size();
    if (N == 0) return 0;

    int idx1     = rng(N);
    int idx1_dag = rng(N);
    int idx2     = rng(N);
    int idx2_dag = rng(N);

    if (not config.try_double_erase(idx1, idx1_dag, idx2, idx2_dag)) return 0; //data is not modified in this case

    int n_fops   = all_d_ops.size();
    return std::pow(N / (n_fops * params.tau_max), 4);
  }

} // namespace inchworm::moves
