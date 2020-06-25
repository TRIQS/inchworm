#include "./solver_core.hpp"

#include "./post_process.hpp"

#include "./mc/measures/u_frame.hpp"
#include "./mc/measures/g_frame.hpp"

#include "./mc/moves/insert.hpp"
#include "./mc/moves/remove.hpp"

#include <triqs/utility/callbacks.hpp>
#include <triqs/utility/macros.hpp>
#include <triqs/mc_tools/mc_generic.hpp>

#define N_STEP 10

namespace inchworm {

  //------------------------------
  // Constructor:
  solver_core::solver_core(constr_params_t const &p) : gf_struct(p.gf_struct), constr_params(p) {

    // Initialize the hybridization function
    Delta_tau = h_tau_t{{p.beta, Fermion, p.n_tau}, p.gf_struct};

    // Determine basis of operators to use
    int n_fops = 0;
    for (auto const &bl : gf_struct) {
      for (auto const &a : bl.second) {
        fops.insert(bl.first, a);
        n_fops++;
      }
    }

    // Setup the linear index map (link Green function structure to fundamental operator set):
    int block_index = 0;
    for (auto const &bl : gf_struct) {
      int inner_index = 0;
      for (auto const &a : bl.second) {
        map_lin_idx_to_block_inner[fops[{bl.first, a}]] = std::make_pair(block_index, inner_index);
        inner_index++;
      }
      block_index++;
    }
  }

  // -------------------------------------------------------------------------------

  //------------------------------
  // FIXME: to be done. For now there is no vanilla solve()
  void solver_core::solve(solve_params_t const &solve_params) {}

  // Common initilization to any solving scheme:
  void solver_core::init(solve_params_t const &solve_params) {

    // Reset the results
    container_set::operator=(container_set{});

    // http://patorjk.com/software/taag/#p=testall&f=Calvin%20S&t=TRIQS%20inchworm%0A
    if (world.rank() == 0)
      std::cout << "\n"
                   "╔╦╗╦═╗╦╔═╗ ╔═╗  ┬┌┐┌┌─┐┬ ┬┬ ┬┌─┐┬─┐┌┬┐\n"
                   " ║ ╠╦╝║║═╬╗╚═╗  │││││  ├─┤││││ │├┬┘│││\n"
                   " ╩ ╩╚═╩╚═╝╚╚═╝  ┴┘└┘└─┘┴ ┴└┴┘└─┘┴└─┴ ┴\n";

    //
    if (solve_params.partition_method != "quantum_numbers")
      TRIQS_RUNTIME_ERROR << "Please use total number for quantum number and use quantum numbers methods for partition of atom_diag";
    //
    h_imp  = solve_params.h_imp;
    ad_imp = {h_imp, fops, solve_params.quantum_numbers};
    u_tau  = make_propagator(ad_imp, constr_params.beta, N_STEP + 1);
    //print_eigensystems(ad_imp);
  }

  //------------------------------
  // solve cthyb (no split point + bare propagator):
  single_step_results_t solver_core::solve_cthyb(solve_params_t const &solve_params, double tau_max) {

    // Merge constr_params and solve_params
    last_solve_params = solve_params;

    // Initialize:
    init(solve_params);

    // use_bare_propagator == cthyb (no split point).
    double tau_split         = 0.0;
    bool use_bare_propagator = true;

    // u_frame recipient:
    u_frame_bare = make_bare_propagator_frame(ad_imp, tau_max, false);
    auto res     = single_step(solve_params, tau_split, tau_max, use_bare_propagator, 0);

    // Finding the ratio between theoretical zeroth order and Monte Carlo sampled zeroth order.
    // Usually first value of u_frame is the most significant, due to order of eigenvalues.
    scalar_t normalization_cte = (double)res.frame_0th_order[0](0, 0) / ((double)u_frame_bare[0](0, 0));
    std::printf("\n\n##################\ncthyb U(tau_max):\n");

    res.normalize(normalization_cte);
    res.print();

    return res;
  } // namespace inchworm

  //------------------------------
  // self consistent solution (one step, with precalculated U(beta) from ED)
  single_step_results_t solver_core::solve_self_consistently(solve_params_t const &solve_params, u_tau_t const &u_tau_, double tau_split,
                                                             double tau_max) {
    // Merge constr_params and solve_params
    last_solve_params = solve_params;

    // Initialize:
    init(solve_params);

    // precalculated propagator:
    u_tau = u_tau_;

    // u_frame recipient:
    u_frame_bare = make_bare_propagator_frame(ad_imp, tau_max, false);

    // calculation of the solution:
    auto res = single_step(solve_params, tau_split, tau_max, false, 0);

    // determine u_frame at zerothr order exactly:
    auto u_frame_zeroth_order = u_tau[0](tau_max - tau_split) * u_tau[0](tau_split);

    // finding the ratio between theoretical zeroth order and Monte Carlo sampled zeroth order:
    scalar_t normalization_cte =
       (double)res.frame_0th_order[0](0, 0) / ((double)u_frame_zeroth_order(0, 0)); //FIXME: need to do better at some point

    std::printf("\n##################\ninchworm U(beta):\n");
    res.normalize(normalization_cte);
    res.print();

    return res;
  } // namespace inchworm

  //------------------------------
  // The inching solution:
  void solver_core::solve_inchworm(solve_params_t const &solve_params) {

    // Initialize:
    init(solve_params);
    beta = constr_params.beta;

    // loop on different inchworm steps
    for (
       int n = 0; n < N_STEP;
       n++) { // FIXME: create a parameter. (at first it was the paramter n_tau, but it is important it is a different one). Now it is just a preprocessor variable. To be done.
      std::printf("\n\ninchworm step %d\n", n);

      // define the tau_split and tau_max for this specific inchworm step.
      //
      // |===================================|----------------------|__________________|
      // 0                                tau_split              tau_max              beta
      //
      // tau_split < tau_max <= beta
      //
      double tau_split = beta * (double)n / (double)N_STEP;
      double tau_max   = beta * (double)(n + 1) / (double)N_STEP;

      // use bare propagator (cthyb) only on the first inchworm iteration:
      bool use_bare_propagator = (n == 0);

      // u_frame recipient:
      u_frame_bare = make_bare_propagator_frame(ad_imp, tau_max, false);

      // calculation of the Monte Carlo solution:
      auto res = single_step(solve_params, tau_split, tau_max, use_bare_propagator, 0);

      // determination of normalization constant:
      scalar_t normalization_cte;
      if (use_bare_propagator) {
        normalization_cte = (double)res.frame_0th_order[0](0, 0) / ((double)u_frame_bare[0](0, 0));
        std::printf("\n\n##################\ncthyb U(tau_max):\n");
      } else {
        auto u_frame_zeroth_order = u_tau[0](tau_max - tau_split) * u_tau[0](tau_split);
        normalization_cte         = (double)res.frame_0th_order[0](0, 0) / ((double)u_frame_zeroth_order(0, 0)); //need to do better at some point
        std::printf("\n\n##################\ninchworm U(tau_max):\n");
      }

      res.normalize(normalization_cte);
      res.print();

      // assign u_frame to the propagator u_tau, in order to be able to use it in next iteration
      assign_frame_to_propagator(u_tau, res.frame, n + 1, 1.);
    }
  } // namespace inchworm

  //------------------------------
  // The Green sampling:
  void solver_core::solve_green(solve_params_t const &solve_params, u_tau_t const &u_tau) {

    // Initialize:
    //exit(1);
    init(solve_params); //FIXME
    beta = constr_params.beta;

    // loop on different inchworm steps
    for (
       int n = 0; n < constr_params.n_tau;
       n++) { // FIXME: create a parameter. (at first it was the paramter n_tau, but it is important it is a different one). Now it is just a preprocessor variable. To be done.
      std::printf("\n\ngreen sampling %d\n", n);

      // define the tau_split and tau_max for this specific inchworm step.
      //
      // d_dag-================================-d-======================|
      // 0                                  tau_split                beta
      //
      // 0 < tau_split < beta
      //
      double tau_split = beta * (double)n / (double)(constr_params.n_tau-1);

      // g_frame recipient:
      // g_frame_bare = make_bare_propagator_frame(ad_imp, tau_max, false);

      auto l = make_u_partial(make_bare_propagator_frame(ad_imp, beta - tau_split, false));
      auto r = make_u_partial(make_bare_propagator_frame(ad_imp, tau_split, false));

      frame_t g_frame_zeroth_order = make_frame(constr_params.gf_struct);

      int n_ops = ad_imp.get_fops().data().size();
      for (int i = 0; i < n_ops; ++i) {
        auto [g_bl, in] = map_lin_idx_to_block_inner.at(i);
        for (int j = 0; j < n_ops; ++j) {
          auto [g_bl_dag, in_dag] = map_lin_idx_to_block_inner.at(j);
          if (g_bl != g_bl_dag) continue;

          // r * ddag_j[bl_idx1](0)
          u_partial_t rddag = apply_op_from_right(r, j, true, ad_imp);

          // l * d_i[bl_idx2](tau)
          u_partial_t ld = apply_op_from_right(l, i, false, ad_imp);

          auto prod = make_u_frame(ld * rddag);

          for (int bl0 = 0; bl0 < ad_imp.n_subspaces(); ++bl0) {
            g_frame_zeroth_order[g_bl](in, in_dag) += trace(prod[bl0]); //FIXME check the order of in and in_dag to be sure.
          }
        }
      }

      TRIQS_PRINT(n); 
      if(n == 0 or n == constr_params.n_tau - 1){
        assign_frame_to_propagator(G_tau, g_frame_zeroth_order, n, 1.);
	continue;
      }

      // calculation of the Monte Carlo solution:
      auto res = single_step(solve_params, tau_split, beta, false, 1);

      // determination of normalization constant:
      scalar_t normalization_cte;

      // FIXME function:
      //--->auto g_frame_zeroth_order = Trace u_tau[0](beta - tau_split) * d_b *  u_tau[0](tau_split) * d_dag_a;

      // FIXME
      normalization_cte = (double)res.frame_0th_order[0](0, 0) / ((double)g_frame_zeroth_order[0](0, 0)); //need to do better at some point
      std::printf("\n\n##################\ninchworm G(tau_split):\n");

      print(res.frame_0th_order);
      getchar();
      res.normalize(normalization_cte);
      res.print();

      assign_frame_to_propagator(G_tau, res.frame, n, 1.);
    }
  } // namespace inchworm

  //------------------------------
  // one Monte Carlo step calculation (common to all solve scheme above):
  single_step_results_t solver_core::single_step(solve_params_t const &solve_params, double tau_split, double tau_max, bool use_bare_propagator,
                                                 int mode) {
    // mode 0 = propagator (inchworm)
    // mode 1 = green function

    params_t params(constr_params, solve_params);
    // Construct the generic Monte-Carlo solver
    triqs::mc_tools::mc_generic<scalar_t> mc(params.random_name, params.random_seed, params.verbosity);

    // Capture random number generator
    auto &rng = mc.get_rng();

    // necessary:
    u_tau_t *u_tau_p = &u_tau;
    if (use_bare_propagator) u_tau_p = nullptr;

    // Create Monte-Carlo configuration
    qmc_config_data_t qmc_config_data{params.gf_struct};

    // Create Monte-Carlo params
    qmc_params_t qmc_params{Delta_tau, map_lin_idx_to_block_inner, ad_imp, u_tau, tau_max, tau_split, use_bare_propagator, mode};

    mc.add_move(moves::insert{qmc_config_data, params, qmc_params, rng}, "insert move");
    mc.add_move(moves::remove{qmc_config_data, params, qmc_params, rng}, "remove move");

    std::vector<long> shape_of_frame;
    if (mode == 0) {
      for (int bl = 0; bl < ad_imp.n_subspaces(); bl++) { shape_of_frame.push_back(ad_imp.get_subspace_dim(bl)); }
    } else {
      for (auto const &[bl, idxlst] : gf_struct) { shape_of_frame.push_back(idxlst.size()); }
    }

    // initialize result container:
    single_step_results_t results{shape_of_frame};
    if (mode == 0)
      results.frame = make_zero_propagator_frame(ad_imp);
    else if (mode == 1)
      results.frame = make_frame(params.gf_struct);
    results.frame_0th_order = results.frame;

    // Register all measurements
    if (mode == 0)
      mc.add_measure(measures::u_frame{params, qmc_config_data, results}, "propagator measurement");
    else if (mode == 1)
      mc.add_measure(measures::g_frame{params, qmc_config_data, results}, "propagator measurement");

    // Perform QMC run and collect results
    mc.warmup_and_accumulate(params.n_warmup_cycles, params.n_cycles, params.length_cycle, triqs::utility::clock_callback(params.max_time));
    mc.collect_results(world);

    // Post Processing
    //if (params.post_process) { post_process(params); }
    return results;
  }

  // -------------------------------------------------------------------------------

  void solver_core::post_process(params_t const &) {

    if (world.rank() == 0)
      std::cout << "\n"
                   "Post-processing ... \n";

    // TODO
  }

  void h5_write(triqs::h5::group h5group, std::string subgroup_name, solver_core const &s) {
    auto grp = h5group.create_group(subgroup_name);
    h5_write_attribute(grp, "TRIQS_HDF5_data_scheme", solver_core::hdf5_scheme());
    h5_write_attribute(grp, "TRIQS_GIT_HASH", std::string(AS_STRING(TRIQS_GIT_HASH)));
    h5_write_attribute(grp, "INCHWORM_GIT_HASH", std::string(AS_STRING(INCHWORM_GIT_HASH)));
    h5_write(grp, "", s.result_set());
    h5_write(grp, "constr_params", s.constr_params);
    h5_write(grp, "last_solve_params", s.last_solve_params);
    //h5_write(grp, "G0_iw", s.G0_iw);
  }

  solver_core solver_core::h5_read_construct(triqs::h5::group h5group, std::string subgroup_name) {
    auto grp           = h5group.open_group(subgroup_name);
    auto constr_params = h5_read<constr_params_t>(grp, "constr_params");
    auto s             = solver_core{constr_params};
    h5_read(grp, "", s.result_set());
    h5_read(grp, "last_solve_params", s.last_solve_params);
    //h5_read(grp, "G0_iw", s.G0_iw);
    return s;
  }

} // namespace inchworm
