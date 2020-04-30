#include "./solver_core.hpp"

#include "./post_process.hpp"

//#include "./mc/measures/sign.hpp"
#include "./mc/measures/u_frame.hpp"
#include "./mc/measures/average_k.hpp"

#include "./mc/moves/insert.hpp"
#include "./mc/moves/remove.hpp"

#include <triqs/utility/callbacks.hpp>
#include <triqs/mc_tools/mc_generic.hpp>

namespace inchworm {

  solver_core::solver_core(constr_params_t const &p) : gf_struct(p.gf_struct), constr_params(p) {

    // Initialize the non-interacting Green function
    //G0_iw = block_gf<imfreq>{{p.beta, Fermion, p.n_iw}, p.gf_struct};

    // Initialize the result containers
    Delta_tau = h_tau_t{{p.beta, Fermion, p.n_tau}, p.gf_struct};
    //G_iw  = G0_iw;

    // determine basis of operators to use
    //fundamental_operator_set fops;
    int n_fops = 0;
    for (auto const &bl : gf_struct) {
      for (auto const &a : bl.second) {
        fops.insert(bl.first, a);
        n_fops++;
      }
    }

    // setup the linear index map
    //std::map<std::pair<int, int>, int> linindex;
    int block_index = 0;
    for (auto const &bl : gf_struct) {
      int inner_index = 0;
      for (auto const &a : bl.second) {
        //linindex[std::make_pair(block_index, inner_index)] = fops[{bl.first, a}];
        map_lin_idx_to_block_inner[fops[{bl.first, a}]] = std::make_pair(block_index, inner_index);
        inner_index++;
      }
      block_index++;
    }

    //Sigma_iw = G0_iw;
  }

  // -------------------------------------------------------------------------------

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
    _h_loc = solve_params.h_int;
    h_diag = {_h_loc, fops, solve_params.quantum_numbers};
    u_tau  = make_propagator(h_diag, constr_params.beta, constr_params.n_tau);
    //print_eigensystems(h_diag);
  }

  void solver_core::solve(solve_params_t const &solve_params) {

    // Merge constr_params and solve_params
    last_solve_params = solve_params;
    init(solve_params);
    //for
    single_step(solve_params, constr_params.beta / 2, constr_params.beta, true);
  }

  void solver_core::solve_single_step(solve_params_t const &solve_params) {

    // Merge constr_params and solve_params
    last_solve_params = solve_params;
    init(solve_params);
    double tau_max   = constr_params.beta;
    double tau_split = constr_params.beta / 2;

    u_frame_bare = make_bare_propagator_frame(h_diag, tau_max, false);
    auto res     = single_step(solve_params, tau_split, tau_max, true);

    double normalization_cte = (double)res.u_frame_0th_order[0](0, 0) / ((double)u_frame_bare[0](0, 0)); //need to do better at some point
    //std::cout << "\n\nnormalization_cte: " << (double)res.u_frame_0th_order[0](0, 0) << "  " << ((double)u_frame_bare[0](0, 0)) << "  " << normalization_cte << "\n";
    std::printf("\n ");
    //for (auto const &B : res.u_frame_0th_order) std::cout << (double) (B/normalization_cte);
    //for (auto &B : res.u_frame_0th_order) {
    //  B /= normalization_cte;
    //  std::cout << B;
    //}
    //std::printf("\n ");
    std::printf("\ncthyb U(beta):\n");
    for (auto &B : res.u_frame) {
      B /= normalization_cte;
      std::cout << B;
    }
    std::cout << "\n\norder: " << res.average_k << "\n";
    for (auto o : res.samples_expansion_order) std::printf("%16d ", o);
    std::printf("\n");
    for (auto o : res.u_expansion_order) std::printf("% 16.5f ", o / normalization_cte);
    std::printf("\n");
  } // namespace inchworm

  void solver_core::solve_self_consistently(solve_params_t const &solve_params, u_tau_t const &u_tau_, double tau_split, double tau_max) {

    // Merge constr_params and solve_params
    last_solve_params = solve_params;
    init(solve_params);
    //double tau_max   = constr_params.beta;
    //double tau_split = constr_params.beta / 2;

    u_tau        = u_tau_;
    u_frame_bare = make_bare_propagator_frame(h_diag, tau_max, false);
    auto res     = single_step(solve_params, tau_split, tau_max, false);

    auto u_frame_zeroth_order = u_tau[0](tau_max - tau_split) * u_tau[0](tau_split);
    double normalization_cte  = (double)res.u_frame_0th_order[0](0, 0) / ((double)u_frame_zeroth_order(0, 0)); //need to do better at some point
    std::printf("\n ");

    std::printf("\ninchworm U_0(beta):\n");
    for (int i = 0; i < res.u_frame.size(); i++) {
      printf("% 4.12f *", (double)u_tau[i](tau_split)(0, 0));
      printf("% 4.12f =", (double)u_tau[i](tau_max - tau_split)(0, 0));
      auto u_frame_zeroth_order = u_tau[i](tau_split) * u_tau[i](tau_max - tau_split);
      printf("% 4.12f\n", (double)(u_frame_zeroth_order(0, 0)));
    }

    std::printf("\n##################\ninchworm U(beta):\n");
    for (int i = 0; i < res.u_frame.size(); i++) { std::cout << (matrix_t)(res.u_frame[i] / normalization_cte); }
    std::cout << "\n\norder: " << res.average_k << "\n";
    for (auto o : res.samples_expansion_order) std::printf("%16d ", o);
    std::printf("\n");
    for (auto o : res.u_expansion_order) std::printf("% 16.5f ", o / normalization_cte);
    std::printf("\n");
  } // namespace inchworm

  void solver_core::solve_inchworm(solve_params_t const &solve_params) {

    // Merge constr_params and solve_params
    //last_solve_params = solve_params;
    init(solve_params);
    //double tau_max   = constr_params.beta;
    //double tau_split = constr_params.beta / 2;
    double beta = constr_params.beta;
    int n_tau   = constr_params.beta;

    for (int n = 0; n < n_tau; n++) {

      double tau_split = beta * n / (n_tau - 1);
      double tau_max   = beta * (n + 1) / (n_tau - 1);

      bool use_bare_propagator = (n == 0 ? true : false);
      u_frame_bare             = make_bare_propagator_frame(h_diag, tau_max, false);
      auto res                 = single_step(solve_params, tau_split, tau_max, use_bare_propagator);

      auto u_frame_zeroth_order = u_tau[0](tau_max - tau_split) * u_tau[0](tau_split);
      double normalization_cte  = (double)res.u_frame_0th_order[0](0, 0) / ((double)u_frame_zeroth_order(0, 0)); //need to do better at some point
      std::printf("\n ");

      std::printf("\ninchworm U(beta):\n");
      for (int i = 0; i < res.u_frame.size(); i++) { std::cout << (matrix_t)(res.u_frame[i] / normalization_cte); }
      std::cout << "\n\norder: " << res.average_k << "\n";
      for (auto o : res.samples_expansion_order) std::printf("%16d ", o);
      std::printf("\n");
      for (auto o : res.u_expansion_order) std::printf("% 16.5f ", o / normalization_cte);
      std::printf("\n");
    }
  } // namespace inchworm

  //------------------------------
  single_step_results_t solver_core::single_step(solve_params_t const &solve_params, double tau_split, double tau_max, bool use_bare_propagator) {

    params_t params(constr_params, solve_params);
    // Construct the generic Monte-Carlo solver
    triqs::mc_tools::mc_generic<scalar_t> mc(params.random_name, params.random_seed, params.verbosity);

    // Capture random number generator
    auto &rng = mc.get_rng();

    // Create Monte-Carlo configuration
    qmc_config_data_t qmc_config_data{h_diag, tau_max, tau_split, &u_tau};

    // Create Monte-Carlo params
    qmc_params_t qmc_params{Delta_tau, map_lin_idx_to_block_inner, h_diag, u_tau, tau_max, tau_split, use_bare_propagator};
    //params, h_diag, u_tau, _Delta_tau, map_lin_idx_to_block_inner};

    mc.add_move(moves::insert{qmc_config_data, qmc_params, rng}, "insert move");
    mc.add_move(moves::remove{qmc_config_data, qmc_params, rng}, "remove move");

    single_step_results_t results(h_diag);
    // Register all measurements
    //mc.add_measure(measures::sign{params, qmc_config_data, results}, "sign measurement");
    mc.add_measure(measures::u_frame{params, qmc_config_data, results}, "propagator measurement");
    mc.add_measure(measures::average_k{params, qmc_config_data, results}, "average perturbation order");

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
