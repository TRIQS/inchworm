#include "./solver_core.hpp"
#include "./u_frame.hpp"
#include "./atom_diag.hpp"
#include "./post_process.hpp"
#include "./measures.hpp"
#include "./moves/insert.hpp"
#include "./moves/remove.hpp"

#include <triqs/utility/callbacks.hpp>
#include <triqs/mc_tools/mc_generic.hpp>

#include <fmt/core.h>

namespace inchworm {

  //------------------------------
  // Constructor:
  solver_core::solver_core(constr_params_t const &cp) : constr_params(cp) {

    // Initialize the hybridization function
    Delta_tau = h_tau_t{{cp.beta, Fermion, cp.n_tau}, cp.gf_struct};

    // Determine basis of operators to use
    fops = fundamental_operator_set{cp.gf_struct};

    // Initialize containers containing all d and d_dag operators by block
    long n_bl = cp.gf_struct.size();
    all_d_ops.resize(n_bl, {});
    all_d_dag_ops.resize(n_bl, {});

    for (auto bl : range(n_bl)) {
      auto [bl_name, bl_size] = cp.gf_struct[bl];
      for (auto idx : range(bl_size)) {
        all_d_ops[bl].push_back({0.0, false, fops[{bl_name, idx}], bl, idx});
        all_d_dag_ops[bl].push_back({0.0, true, fops[{bl_name, idx}], bl, idx});
      }
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

    if (sp.partition_method == "automatic") {
      ASSERT(sp.quantum_numbers.empty());
      ad_imp = {sp.h_imp, create_effective_hyb(constr_params.gf_struct), fops}; // FIXME Change order of arguments?
    } else if (sp.partition_method == "quantum_numbers") {
      ad_imp = {sp.h_imp, fops, sp.quantum_numbers};
    } else {
      TRIQS_RUNTIME_ERROR
         << "Unknown partition method! Please choose 'automatic' or 'quantum_number' and set solve_params.quantum_numbers accordingly";
    }
  }

  //------------------------------
  // solve cthyb (no split point + bare propagator):
  qmc_results_t solver_core::solve_cthyb(solve_params_t const &solve_params, double tau_max) {

    // Initialize solver
    this->init(solve_params);

    // use_bare_propagator == cthyb (no split point).
    double tau_split         = 0.0;
    bool use_bare_propagator = true;

    // Execute cthyb sampling
    auto res = qmc_step(solve_params, tau_split, tau_max, use_bare_propagator, MODE::PROPAGATOR);

    // Finding the ratio between theoretical zeroth order and Monte Carlo sampled zeroth order.
    // Usually first value of u_frame is the most significant, due to order of eigenvalues.
    auto u_frame_bare          = make_bare_u_frame(ad_imp, tau_max);
    scalar_t normalization_cte = frobenius_norm(res.frame_0th_order) / frobenius_norm(u_frame_bare);
    if (normalization_cte == 0)
      TRIQS_RUNTIME_ERROR << "Failed to calculate normalization ratio due to insufficient sampling of zeroth order propagator";
    res.normalize(normalization_cte);

    // Print results
    res.print(solve_params.verbosity);

    return res;
  }

  //------------------------------

  // Self consistent solution (one step, with precalculated U(beta) from ED)
  qmc_results_t solver_core::solve_self_consistently(solve_params_t const &solve_params, u_tau_t const &u_tau_, double tau_split, double tau_max) {
    // Initialize solver
    this->init(solve_params);

    // precalculated propagator:
    u_tau = u_tau_;

    // Initialize the interpolator
    u_interpolator = interpolator_t(u_tau, u_tau[0].mesh().size());

    // Execute u_tau self-consistency sampling
    auto res = qmc_step(solve_params, tau_split, tau_max, false, MODE::PROPAGATOR);

    // Normalize the result using the ratio between theoretical zeroth order and Monte Carlo sampled zeroth order
    auto u_frame_zeroth_order  = frame_t{u_interpolator(tau_max - tau_split) * u_interpolator(tau_split)};
    scalar_t normalization_cte = frobenius_norm(res.frame_0th_order) / frobenius_norm(u_frame_zeroth_order);
    if (normalization_cte == 0)
      TRIQS_RUNTIME_ERROR << "Failed to calculate normalization ratio due to insufficient sampling of zeroth order propagator";
    res.normalize(normalization_cte);

    // Print results
    res.print(solve_params.verbosity);

    return res;
  }

  //------------------------------

  void solver_core::solve_inchworm(solve_params_t const &solve_params, bool use_cthyb) {

    // Initialize solver
    this->init(solve_params);

    // Reset the results
    container_set::operator=(container_set{});

    // Initialize empty propagator
    auto u_tau_zero = u_tau_t{{constr_params.beta, Fermion, constr_params.n_tau_inch}, ad_imp.get_subspace_dims()};
    u_tau_zero()    = 0.;
    u_tau           = u_tau_zero;
    for (auto &ubl : u_tau) {
      for (int i = 0; i < ubl.target_shape()[0]; ++i) ubl[0](i, i) = 1;
    }
    u_tau_by_order = {u_tau};

    if (solve_params.verbosity > 0) std::cout << "\nStarting inchworm calculation of the propagator.. \n";

    // loop on different inchworm steps
    for (auto n : range(1, constr_params.n_tau_inch)) {
      if (solve_params.verbosity > 0) std::printf("\n ..step %ld/%d\n", n, constr_params.n_tau_inch - 1);

      if (solve_params.n_tau_inch_stop && n > *solve_params.n_tau_inch_stop) break;

      // define the tau_split and tau_max for this specific inchworm step.
      //
      // |===================================|----------------------|__________________|
      // 0                                tau_split              tau_max              beta
      //
      // tau_split < tau_max <= beta
      //
      double tau_split = use_cthyb ? 0.0 : u_tau[0].mesh().index_to_point(n - 1);
      double tau_max   = u_tau[0].mesh().index_to_point(n);

      // use bare propagator (cthyb) only on the first inchworm iteration:
      bool use_bare_propagator = (n == 1) or use_cthyb;

      // Initialize the interpolator
      if (not use_bare_propagator) { u_interpolator = interpolator_t(u_tau, n); }

      // calculation of the Monte Carlo solution:
      auto res = qmc_step(solve_params, tau_split, tau_max, use_bare_propagator, MODE::PROPAGATOR);

      // Normalize the result using the ratio between theoretical zeroth order and Monte Carlo sampled zeroth order
      scalar_t normalization_cte;
      if (use_bare_propagator) {
        auto u_frame_bare = make_bare_u_frame(ad_imp, tau_max);
        normalization_cte = frobenius_norm(res.frame_0th_order) / frobenius_norm(u_frame_bare);
      } else {
        auto u_frame_zeroth_order = frame_t{u_interpolator(tau_max - tau_split) * u_interpolator(tau_split)};
        normalization_cte         = frobenius_norm(res.frame_0th_order) / frobenius_norm(u_frame_zeroth_order);
      }
      if (normalization_cte == 0)
        TRIQS_RUNTIME_ERROR << "Failed to calculate normalization ratio due to insufficient sampling of zeroth order propagator";
      res.normalize(normalization_cte);

      // Print results
      res.print(solve_params.verbosity);

      // Assign the current frame to the propagator u_tau (will be used in the next iteration)
      set_frame(res.frame, u_tau, n);

      // Assign u_frame_by_order to the propagator u_tau_by_order
      if (solve_params.measure_frame_by_order) {
        if (res.frame_by_order.size() > u_tau_by_order.size()) u_tau_by_order.resize(res.frame_by_order.size(), u_tau_zero);
        for (auto k : range(res.frame_by_order.size())) set_frame(res.frame_by_order[k], u_tau_by_order[k], n);
      }

      // Initialize other results
      if (solve_params.measure_order_histogram) order_histograms.push_back(res.order_histogram);

      // Break out of loop when interrupted by signal
      if (res.status == 2) break;
    }
  }

  //------------------------------

  void solver_core::solve_green(solve_params_t const &solve_params) {

    // Initialize solver
    this->init(solve_params);

    // Reset the results
    container_set::operator=(container_set{});

    if (solve_params.verbosity > 0) std::cout << "\nStarting Green function calculation.. \n";

    // Initialize the interpolator
    u_interpolator = interpolator_t(u_tau, constr_params.n_tau_inch);

    double beta = constr_params.beta;

    // Initialize the Green function container
    G_tau = g_tau_t{{beta, Fermion, constr_params.n_tau_green}, constr_params.gf_struct};

    // Calculate Tr U(beta)
    scalar_t Tr_Ubeta = 0.0;
    for (int bl = 0; bl < u_tau.size(); bl++) Tr_Ubeta += trace(u_tau[bl](beta));

    // Treat n == 0 and n == n_tau -1 seperately
    frame_t g_frame_n0 = make_bare_g_frame(ad_imp, u_tau, constr_params.gf_struct, 0.0, beta) / Tr_Ubeta;
    frame_t g_frame_nB = make_bare_g_frame(ad_imp, u_tau, constr_params.gf_struct, beta, beta) / Tr_Ubeta;
    set_frame(g_frame_n0, G_tau, 0);
    set_frame(g_frame_nB, G_tau, constr_params.n_tau_green - 1);

    // loop on different green function tau values
    // n == 0 and n == n_tau - 1 already treated
    for (auto n : range(1, constr_params.n_tau_green - 1)) {
      if (solve_params.verbosity > 0) std::printf(" ..step %ld/%d\n", n, constr_params.n_tau_green - 2);

      // define the tau_split and tau_max for this specific inchworm step.
      //
      // d_dag-================================-d-======================|
      // 0                                  tau_split                beta
      //
      // 0 < tau_split < beta
      //
      double tau_split = G_tau[0].mesh().index_to_point(n);

      // calculation of the Monte Carlo solution:
      auto res = qmc_step(solve_params, tau_split, beta, false, MODE::GREENFUNCTION);

      // Normalize the result using the ratio between theoretical zeroth order and Monte Carlo sampled zeroth order
      frame_t g_frame_zeroth_order = make_bare_g_frame(ad_imp, u_tau, constr_params.gf_struct, tau_split, beta);
      scalar_t normalization_cte   = Tr_Ubeta * frobenius_norm(res.frame_0th_order) / frobenius_norm(g_frame_zeroth_order);
      if (normalization_cte == 0)
        TRIQS_RUNTIME_ERROR << "Failed to calculate normalization ratio due to insufficient sampling of zeroth order Green function";
      res.normalize(normalization_cte);

      // Print results
      res.print(solve_params.verbosity);

      // Assign the current frame to G_tau
      set_frame(res.frame, G_tau, n);

      // Initialize other results
      if (solve_params.measure_order_histogram) order_histograms.push_back(res.order_histogram);

      // Break out of loop when interrupted by signal
      if (res.status == 2) break;
    }
  }

  //------------------------------

  // One Monte Carlo sampling step (common to all solve scheme above):
  qmc_results_t solver_core::qmc_step(solve_params_t const &solve_params, double tau_split, double tau_max, bool use_bare_propagator, MODE mode) {
    params_t params(constr_params, solve_params);

    // Construct the generic Monte-Carlo solver
    triqs::mc_tools::mc_generic<scalar_t> mc(params.random_name, params.random_seed, params.verbosity);

    // Capture random number generator
    auto &rng = mc.get_rng();

    // Create frame and Monte-Carlo configuration
    auto frame =
       (mode == MODE::PROPAGATOR) ? make_bare_u_frame(ad_imp, tau_max) : make_bare_g_frame(ad_imp, u_tau, params.gf_struct, tau_split, tau_max);
    config_t config(frame, params.gf_struct, {0.0, tau_split});

    // Create Monte-Carlo params
    qmc_params_t qmc_params{tau_max, tau_split, use_bare_propagator, mode};

    // Initialize the accumulators for the tau difference statistics
    long n_bl = params.gf_struct.size();
    std::vector<nda::array<accumulator<double>, 3>> tau_diff_stat(n_bl);
    for (auto bl : range(n_bl)) {
      tau_diff_stat[bl] = nda::array<accumulator<double>, 3>{params.gf_struct[bl].second, 2, 2};
      tau_diff_stat[bl] = accumulator<double>{0.0, 0, -1}; // Lin-binning only
    }

    // FIXME Reset operator widths manually
    for (auto bl : range(n_bl)) {
      for (auto &op : all_d_ops[bl]) {
        op.left_width  = fop_t{}.left_width;
        op.right_width = fop_t{}.right_width;
      }
      for (auto &op : all_d_dag_ops[bl]) {
        op.left_width  = fop_t{}.left_width;
        op.right_width = fop_t{}.right_width;
      }
    }

    // Add moves
    using moves::KIND;
    mc.add_move(moves::insert<KIND::Single>{config, frame, qmc_params, *this, rng, tau_diff_stat}, "insert move");
    mc.add_move(moves::remove<KIND::Single>{config, frame, qmc_params, *this, rng}, "remove move");

    if (params.use_double_insertion) {
      mc.add_move(moves::insert<KIND::Double>{config, frame, qmc_params, *this, rng, tau_diff_stat}, "double insert move", 0.5);
      mc.add_move(moves::remove<KIND::Double>{config, frame, qmc_params, *this, rng}, "double remove move", 0.5);
      mc.add_move(moves::insert<KIND::DoubleEqBl>{config, frame, qmc_params, *this, rng, tau_diff_stat}, "double insert move equal blocks", 0.5);
      mc.add_move(moves::remove<KIND::DoubleEqBl>{config, frame, qmc_params, *this, rng}, "double remove move equal blocks", 0.5);
    }

    // Initialize result container
    std::vector<int> shape_of_frame;
    if (mode == MODE::PROPAGATOR) {
      for (int bl = 0; bl < ad_imp.n_subspaces(); bl++) { shape_of_frame.push_back(ad_imp.get_subspace_dim(bl)); }
    } else { // MODE::GREENFUNCTION
      for (auto const &[blname, blsize] : params.gf_struct) { shape_of_frame.push_back(blsize); }
    }
    auto results = qmc_results_t{shape_of_frame};

    // ----- Run the warmup
    if (params.verbosity > 0) std::printf("     Warming up ...\n");
    moves::base_move::reweighting_cutoff = 0;
    moves::base_move::reweighting_coeffs.clear();
    auto length_cycle = params.length_cycle.value_or(1);
    int status        = mc.warmup(params.n_warmup_cycles, length_cycle, triqs::utility::clock_callback(params.max_time));

    // ----- Define Function for measuring and setting the operator distribution widths
    bool tau_diff_stat_done = false;
    auto gather_tau_diff_stat = [&]() {
      // Gather the tau-diff statistics
      moves::base_move::gather_tau_diff_stat = true;
      int status                          = mc.warmup(params.n_warmup_cycles, length_cycle, triqs::utility::clock_callback(params.max_time));
      moves::base_move::gather_tau_diff_stat = false;

      // Adjust the operators accordingly
      long below_threshold_count = 0;
      for (auto bl : range(n_bl)) {
        for (auto &op : all_d_ops[bl]) { // FIXME join(all_d_ops[bl], all_d_dag_ops[bl])
          if (mpi::all_reduce(tau_diff_stat[bl](op.idx, op.dag, 0).n_lin_bins()) > 0)
            op.left_width = mean_mpi(world, tau_diff_stat[bl](op.idx, op.dag, 0).linear_bins());
          else
            ++below_threshold_count;

          if (mpi::all_reduce(tau_diff_stat[bl](op.idx, op.dag, 1).n_lin_bins()) > 0)
            op.right_width = mean_mpi(world, tau_diff_stat[bl](op.idx, op.dag, 1).linear_bins());
          else
            ++below_threshold_count;

          if (params.verbosity > 0) PRINT(op);
        }
        for (auto &op : all_d_dag_ops[bl]) {
          if (mpi::all_reduce(tau_diff_stat[bl](op.idx, op.dag, 0).n_lin_bins()) > 0)
            op.left_width = mean_mpi(world, tau_diff_stat[bl](op.idx, op.dag, 0).linear_bins());
          else
            ++below_threshold_count;

          if (mpi::all_reduce(tau_diff_stat[bl](op.idx, op.dag, 1).n_lin_bins()) > 0)
            op.right_width = mean_mpi(world, tau_diff_stat[bl](op.idx, op.dag, 1).linear_bins());
          else
            ++below_threshold_count;

          if (params.verbosity > 0) PRINT(op);
        }
        tau_diff_stat_done = true;
      }
      if (below_threshold_count > 0 && params.verbosity > 0)
        fmt::print("         Found {} operators with insertion count below threshold, please raise n_warmup_cycles\n", below_threshold_count);
    };

    // -----  Run the callibration loop
    size_t hist_max_idx = 0;
    if (params.verbosity > 0) {
      std::printf("     Callibrating ...\n");
      std::printf("         %-12s| %-12s| %-12s| %-12s| %-12s| %-16s\n", "hist0", "autocorr", "acc insert", "acc remove", "new coeff0",
                  "new length_cycle");
    }
    for (int n = 1; status == 0; ++n) {
      if (params.verbosity > 2) std::cout << "\nCallibration-loop " << n << "\n";

      auto callibration_results = results;
      mc.add_measure(measures::average_order{params, config, callibration_results}, "measure the average perturbation order");
      mc.add_measure(measures::order_histogram{params, config, callibration_results}, "measure the perturbation order histogram");
      mc.add_measure(measures::autocorr{params, config, frame, callibration_results}, "measure the autocorrelation time");

      status = mc.accumulate(params.n_callibration_cycles, length_cycle, triqs::utility::clock_callback(params.max_time));
      if (status != 0) break;
      mc.collect_results(world);
      mc.clear_measures();

      // Update reweighting coefficients based on the perturbation order histogram taking the histogram max as a reference
      auto const &hist = callibration_results.order_histogram;
      if (n == 1) { // Set the hist_max_idx and reweighting cutoff on the first iteration
        hist_max_idx                         = std::distance(begin(hist), max_element(begin(hist), end(hist)));
        moves::base_move::reweighting_cutoff = std::max(1ul, hist_max_idx);
        moves::base_move::reweighting_coeffs.resize(moves::base_move::reweighting_cutoff, 1.0);
      }
      for (auto k : range(hist_max_idx)) { // Scale up the weight for orders below hist_max_idx
        auto hist_ratio = std::max(1.0, hist[hist_max_idx] / std::max(hist[k], 0.5 / params.n_callibration_cycles));
        moves::base_move::reweighting_coeffs[k] *= hist_ratio;
        if (config.size() == k) config.imp_weight *= hist_ratio;
      }
      if (params.max_prob_zeroth_order < hist[0]) { // Scale down the weight for the zeroth order if necessary
        auto zero_reweight = 0.95 * std::max(1.0 - hist[0], 0.5 / params.n_callibration_cycles) / hist[0] * params.max_prob_zeroth_order
           / (1.0 - params.max_prob_zeroth_order);
        moves::base_move::reweighting_coeffs[0] *= zero_reweight;
        if (config.size() == 0) config.imp_weight *= zero_reweight;
      }

      // When not reweighting, auto-deduce cycle length if not set
      if (hist[0] > 0.7 * hist[hist_max_idx] and hist[0] <= params.max_prob_zeroth_order)
        length_cycle = params.length_cycle.value_or(1.0 + std::ceil(length_cycle * std::pow(0.1 + callibration_results.auto_corr_time, 0.85)));

      if (params.verbosity > 0) {
        auto acc_rates = mc.get_acceptance_rates();
        std::printf("         %-12.3f| %-12.3f| %-12.3f| %-12.3f| %-12.3e| %-16d\n", hist[0], callibration_results.auto_corr_time,
                    acc_rates.at("insert move"), acc_rates.at("remove move"), moves::base_move::reweighting_coeffs[0], length_cycle);
      }

      // When not reweighting, set operator tau distribution widths
      if (hist[0] > 0.7 * hist[hist_max_idx] and hist[0] <= params.max_prob_zeroth_order and not tau_diff_stat_done) {
	gather_tau_diff_stat();
	continue;
      }

      // Iterate the callibration until the zeroth order is sampled with finite probability
      if (hist[0] > 0.7 * hist[hist_max_idx] and hist[0] <= params.max_prob_zeroth_order
          and (params.length_cycle.has_value() or callibration_results.auto_corr_time < 5.0))
        break;
    }

    // Register all measurements
    mc.add_measure(measures::frame{params, config, frame, results}, "measure the propagator / green function frame");
    if (params.measure_average_order) mc.add_measure(measures::average_order{params, config, results}, "measure the average perturbation order");
    if (params.measure_order_histogram)
      mc.add_measure(measures::order_histogram{params, config, results}, "measure the perturbation order histogram");
    if (params.measure_frame_by_order)
      mc.add_measure(measures::frame_by_order{params, config, frame, results}, "measure the propagator / green function frame by order");

    // Perform QMC run and collect results
    if (status == 0) {
      if (params.verbosity > 0) std::printf("     Accumulating ...\n");
      results.status = mc.accumulate(params.n_cycles, length_cycle, triqs::utility::clock_callback(params.max_time));
      mc.collect_results(world);
      if (params.max_order && results.order_histogram[*params.max_order] > 0.0)
        if (world.rank() == 0)
          std::cout << "WARNING: Maximum perturbation order was sampled with a finite probability of " << results.order_histogram[*params.max_order]
                    << ". Check convergence w.r.t. max_order!\n";
    }

    // Post Processing
    //if (params.post_process) { post_process(params); }

    if (results.status == 2) std::cerr << "Warning: Inchworm was interrupted by signal\n";

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
    h5_write(grp, "ad_imp", s.ad_imp);
    h5_write(grp, "Delta_tau", s.Delta_tau);
    h5_write(grp, "u_tau", s.u_tau);
  }

  solver_core solver_core::h5_read_construct(h5::group h5group, std::string subgroup_name) {
    auto grp           = h5group.open_group(subgroup_name);
    auto constr_params = h5_read<constr_params_t>(grp, "constr_params");
    auto s             = solver_core{constr_params};
    h5_read(grp, "", s.result_set());
    h5_read(grp, "last_solve_params", s.last_solve_params);
    h5_read(grp, "ad_imp", s.ad_imp);
    h5_read(grp, "Delta_tau", s.Delta_tau);
    h5_read(grp, "u_tau", s.u_tau);
    return s;
  }

} // namespace inchworm
