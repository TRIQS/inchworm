#include "./solver_core.hpp"

#include "./post_process.hpp"

#include "./mc/measures/u_frame.hpp"
#include "./mc/measures/g_frame.hpp"

#include "./mc/moves/insert.hpp"
#include "./mc/moves/remove.hpp"
#include "./mc/moves/double_insert.hpp"
#include "./mc/moves/double_remove.hpp"

#include <triqs/utility/callbacks.hpp>
#include <triqs/utility/macros.hpp>
#include <triqs/mc_tools/mc_generic.hpp>

namespace inchworm {

  //------------------------------
  // Constructor:
  solver_core::solver_core(constr_params_t const &cp) : constr_params(cp) {

    // Initialize the hybridization function
    Delta_tau = h_tau_t{{cp.beta, Fermion, cp.n_tau}, cp.gf_struct};

    // Determine basis of operators to use
    int n_fops = 0;
    for (auto const &[blname, blsize] : cp.gf_struct) {
      for (auto const &a : range(blsize)) {
        fops.insert(blname, a);
        n_fops++;
      }
    }

    // Setup the linear index map (link Green function structure to fundamental operator set):
    int block_index = 0;
    for (auto const &[blname, blsize] : cp.gf_struct) {
      for (int idx : range(blsize)) { map_lin_idx_to_block_inner[fops[{blname, idx}]] = std::make_pair(block_index, idx); }
      block_index++;
    }
  }

  // -------------------------------------------------------------------------------

  //------------------------------
  // High-level solve function. Run inchworm first to get u_tau
  // and then sample the Green function
  void solver_core::solve(solve_params_t const &solve_params) {

    // http://patorjk.com/software/taag/#p=testall&f=Calvin%20S&t=TRIQS%20inchworm%0A
    if (solve_params.verbosity > 0)
      std::cout << "\n"
                   "╔╦╗╦═╗╦╔═╗ ╔═╗  ┬┌┐┌┌─┐┬ ┬┬ ┬┌─┐┬─┐┌┬┐\n"
                   " ║ ╠╦╝║║═╬╗╚═╗  │││││  ├─┤││││ │├┬┘│││\n"
                   " ╩ ╩╚═╩╚═╝╚╚═╝  ┴┘└┘└─┘┴ ┴└┴┘└─┘┴└─┴ ┴\n";

    // Run inchworm to calculate S.u_tau
    solve_inchworm(solve_params);

    // Sample the Green function S.G_tau
    solve_green(solve_params);
  }

  // Common initilization to any solving scheme:
  void solver_core::init(solve_params_t const &sp) {

    // Store solve_params
    last_solve_params = sp;

    // Reset the results
    container_set::operator=(container_set{});

    // -- Atom Diag Object
    if (sp.partition_method != "quantum_numbers")
      TRIQS_RUNTIME_ERROR << "Please use total number for quantum number and use quantum numbers methods for partition of atom_diag";

    if (sp.quantum_numbers.empty()) {
      ad_imp = {sp.h_imp, create_effective_hyb(constr_params.gf_struct), fops}; // Change order of arguments?
    } else {
      ad_imp = {sp.h_imp, fops, sp.quantum_numbers};
    }
  }

  //------------------------------
  // solve cthyb (no split point + bare propagator):
  single_step_results_t solver_core::solve_cthyb(solve_params_t const &solve_params, double tau_max) {

    // Initialize solver
    this->init(solve_params);

    // use_bare_propagator == cthyb (no split point).
    double tau_split         = 0.0;
    bool use_bare_propagator = true;

    // Execute cthyb sampling
    auto res = single_step(solve_params, tau_split, tau_max, use_bare_propagator, 0);

    // Finding the ratio between theoretical zeroth order and Monte Carlo sampled zeroth order.
    // Usually first value of u_frame is the most significant, due to order of eigenvalues.
    auto u_frame_bare          = make_bare_propagator_frame(ad_imp, tau_max, false);
    scalar_t normalization_cte = (double)res.frame_0th_order[0](0, 0) / ((double)u_frame_bare[0](0, 0));
    res.normalize(normalization_cte);

    // Print results
    if (solve_params.verbosity > 3) {
      std::printf("\n\n##################\ncthyb U(tau_max):\n");
      res.print(solve_params.verbosity);
    }

    return res;
  }

  //------------------------------

  // Self consistent solution (one step, with precalculated U(beta) from ED)
  single_step_results_t solver_core::solve_self_consistently(solve_params_t const &solve_params, u_tau_t const &u_tau_, double tau_split,
                                                             double tau_max) {
    // Initialize solver
    this->init(solve_params);

    // precalculated propagator:
    u_tau = u_tau_;

    // Execute u_tau self-consistency sampling
    auto res = single_step(solve_params, tau_split, tau_max, false, 0);

    // Normalize the result using the ratio between theoretical zeroth order and Monte Carlo sampled zeroth order
    auto u_frame_zeroth_order  = u_tau[0](tau_max - tau_split) * u_tau[0](tau_split);
    scalar_t normalization_cte = (double)res.frame_0th_order[0](0, 0) / ((double)u_frame_zeroth_order(0, 0)); //FIXME: need to do better at some point
    res.normalize(normalization_cte);

    // Print results
    if (solve_params.verbosity > 3) {
      std::printf("\n##################\ninchworm U(beta):\n");
      res.print(solve_params.verbosity);
    }

    return res;
  }

  //------------------------------

  void solver_core::solve_inchworm(solve_params_t const &solve_params) {

    // Initialize solver
    this->init(solve_params);

    // Initialize empty propagator
    u_tau = make_propagator(ad_imp, constr_params.beta, constr_params.n_tau_inch);

    if (solve_params.verbosity > 0) std::cout << "\nStarting inchworm calculation of the propagator.. \n";

    double beta = constr_params.beta;

    int n_step = constr_params.n_tau_inch - 1;
    // loop on different inchworm steps
    for (int n = 0; n < n_step; n++) {
      if (solve_params.verbosity > 0) std::printf(" ..step %d/%d\n", n + 1, n_step);

      // define the tau_split and tau_max for this specific inchworm step.
      //
      // |===================================|----------------------|__________________|
      // 0                                tau_split              tau_max              beta
      //
      // tau_split < tau_max <= beta
      //
      double tau_split = beta * (double)n / (double)n_step;
      double tau_max   = beta * (double)(n + 1) / (double)n_step;

      // use bare propagator (cthyb) only on the first inchworm iteration:
      bool use_bare_propagator = (n == 0);

      // calculation of the Monte Carlo solution:
      auto res = single_step(solve_params, tau_split, tau_max, use_bare_propagator, 0);

      // Normalize the result using the ratio between theoretical zeroth order and Monte Carlo sampled zeroth order
      scalar_t normalization_cte;
      if (use_bare_propagator) {
        auto u_frame_bare = make_bare_propagator_frame(ad_imp, tau_max, false);
        normalization_cte = (double)res.frame_0th_order[0](0, 0) / ((double)u_frame_bare[0](0, 0));
      } else {
        auto u_frame_zeroth_order = u_tau[0](tau_max - tau_split) * u_tau[0](tau_split);
        normalization_cte         = (double)res.frame_0th_order[0](0, 0) / ((double)u_frame_zeroth_order(0, 0)); //need to do better at some point
      }
      res.normalize(normalization_cte);

      // Print results
      if (solve_params.verbosity > 3) {
        std::printf("\n\n##################\ninchworm U(tau_max):\n");
        res.print(solve_params.verbosity);
      }

      // Assign u_frame to the propagator u_tau, in order to be able to use it in next iteration
      assign_frame_to_propagator(u_tau, res.frame, n + 1, 1.);
    }
  }

  //------------------------------

  void solver_core::solve_green(solve_params_t const &solve_params) {

    // Initialize solver
    this->init(solve_params);

    if (solve_params.verbosity > 0) std::cout << "\nStarting Green function calculation.. \n";

    double beta = constr_params.beta;

    // Initialize the Green function container
    G_tau = g_tau_t{{beta, Fermion, constr_params.n_tau_green}, constr_params.gf_struct};

    // --- Treat n == 0 and n == n_tau -1 seperately

    frame_t g_frame_n0 = make_bare_g_frame(ad_imp, u_tau, map_lin_idx_to_block_inner, constr_params.gf_struct, 0.0, beta);
    frame_t g_frame_nB = make_bare_g_frame(ad_imp, u_tau, map_lin_idx_to_block_inner, constr_params.gf_struct, beta, beta);

    // Calculate Tr U(beta)
    scalar_t Tr_Ubeta = 0.0;
    for (int bl = 0; bl < u_tau.size(); bl++) Tr_Ubeta += trace(u_tau[bl](beta));

    assign_frame_to_propagator(G_tau, g_frame_n0, 0, 1. / Tr_Ubeta);
    assign_frame_to_propagator(G_tau, g_frame_nB, constr_params.n_tau_green - 1, 1. / Tr_Ubeta);

    // loop on different inchworm steps
    // n == 0 and n == n_tau - 1 already treated
    for (int n = 1; n < constr_params.n_tau_green - 1; n++) {
      if (solve_params.verbosity > 0) std::printf(" ..step %d/%d\n", n, constr_params.n_tau_green - 2);

      // define the tau_split and tau_max for this specific inchworm step.
      //
      // d_dag-================================-d-======================|
      // 0                                  tau_split                beta
      //
      // 0 < tau_split < beta
      //
      double tau_split = beta * (double)n / (double)(constr_params.n_tau_green - 1);

      // calculation of the Monte Carlo solution:
      auto res = single_step(solve_params, tau_split, beta, false, 1);

      // Normalize the result using the ratio between theoretical zeroth order and Monte Carlo sampled zeroth order
      frame_t g_frame_zeroth_order = make_bare_g_frame(ad_imp, u_tau, map_lin_idx_to_block_inner, constr_params.gf_struct, tau_split, beta);
      scalar_t normalization_cte   = Tr_Ubeta * (double)res.frame_0th_order[0](0, 0) / ((double)g_frame_zeroth_order[0](0, 0));
      res.normalize(normalization_cte);

      // Print results
      if (solve_params.verbosity > 3) {
        std::printf("\n\n##################\ninchworm G(tau_split):\n");
        res.print(solve_params.verbosity);
      }

      assign_frame_to_propagator(G_tau, res.frame, n, 1.);
    }
  }

  //------------------------------

  // One Monte Carlo sampling step (common to all solve scheme above):
  single_step_results_t solver_core::single_step(solve_params_t const &solve_params, double tau_split, double tau_max, bool use_bare_propagator,
                                                 int mode) {
    // mode 0 = propagator (inchworm)
    // mode 1 = green function

    params_t params(constr_params, solve_params);

    // Construct the generic Monte-Carlo solver
    triqs::mc_tools::mc_generic<scalar_t> mc(params.random_name, params.random_seed, params.verbosity);

    // Capture random number generator
    auto &rng = mc.get_rng();

    // Create Monte-Carlo configuration
    qmc_config_data_t qmc_config_data{params.gf_struct};
    if (mode == 0) {
      qmc_config_data.u_partial = make_u_partial(make_bare_propagator_frame(ad_imp, tau_split, false));
    } else {
      qmc_config_data.g_frame = make_bare_g_frame(ad_imp, u_tau, map_lin_idx_to_block_inner, params.gf_struct, tau_split, params.beta);
    }

    // Create Monte-Carlo params
    qmc_params_t qmc_params{Delta_tau, map_lin_idx_to_block_inner, ad_imp, u_tau, tau_max, tau_split, use_bare_propagator, mode};

    // Add moves
    mc.add_move(moves::insert{qmc_config_data, params, qmc_params, rng}, "insert move");
    mc.add_move(moves::remove{qmc_config_data, params, qmc_params, rng}, "remove move");

    mc.add_move(moves::double_insert{qmc_config_data, params, qmc_params, rng}, "double insert move");
    mc.add_move(moves::double_remove{qmc_config_data, params, qmc_params, rng}, "double remove move");

    // Determine the shape of the result
    std::vector<long> shape_of_frame;
    if (mode == 0) {
      for (int bl = 0; bl < ad_imp.n_subspaces(); bl++) { shape_of_frame.push_back(ad_imp.get_subspace_dim(bl)); }
    } else {
      for (auto const &[blname, blsize] : params.gf_struct) { shape_of_frame.push_back(blsize); }
    }

    // Initialize result container
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
      mc.add_measure(measures::g_frame{params, qmc_config_data, results}, "green function measurement");

    // Perform QMC run and collect results
    int status =
       mc.warmup_and_accumulate(params.n_warmup_cycles, params.n_cycles, params.length_cycle, triqs::utility::clock_callback(params.max_time));
    mc.collect_results(world);

    // Post Processing
    //if (params.post_process) { post_process(params); }

    if (status == 2) TRIQS_RUNTIME_ERROR << "Inchworm was stopped by signal";

    return results;
  }

  // -------------------------------------------------------------------------------

  void solver_core::post_process(params_t const &) {

    if (world.rank() == 0)
      std::cout << "\n"
                   "Post-processing ... \n";

    // TODO
  }

  void h5_write(h5::group h5group, std::string subgroup_name, solver_core const &s) {
    auto grp = h5group.create_group(subgroup_name);
    h5_write_attribute(grp, "Format", solver_core::hdf5_format());
    h5_write_attribute(grp, "TRIQS_GIT_HASH", std::string(AS_STRING(TRIQS_GIT_HASH)));
    h5_write_attribute(grp, "INCHWORM_GIT_HASH", std::string(AS_STRING(INCHWORM_GIT_HASH)));
    h5_write(grp, "", s.result_set());
    h5_write(grp, "constr_params", s.constr_params);
    h5_write(grp, "last_solve_params", s.last_solve_params);
    h5_write(grp, "u_tau", s.u_tau);
  }

  solver_core solver_core::h5_read_construct(h5::group h5group, std::string subgroup_name) {
    auto grp           = h5group.open_group(subgroup_name);
    auto constr_params = h5_read<constr_params_t>(grp, "constr_params");
    auto s             = solver_core{constr_params};
    h5_read(grp, "", s.result_set());
    h5_read(grp, "last_solve_params", s.last_solve_params);
    h5_read(grp, "u_tau", s.u_tau);
    return s;
  }

} // namespace inchworm
