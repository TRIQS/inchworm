/*******************************************************************************
 *
 * inchworm: A TRIQS based impurity solver
 *
 * Copyright (c) 2019 The Simons foundation
 *   authors: Nils Wentzell
 *
 * inchworm is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * inchworm is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * inchworm. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#include "./solver_core.hpp"

#include "./post_process.hpp"
#include "./measures/simple.hpp"
#include "./moves/update.hpp"

#include <triqs/utility/callbacks.hpp>
#include <triqs/mc_tools/mc_generic.hpp>

namespace inchworm {

  solver_core::solver_core(constr_params_t const &p) : constr_params(p), gf_struct(p.gf_struct) {

    // Initialize the non-interacting Green function
    G0_iw = block_gf<imfreq>{{p.beta, Fermion, p.n_iw}, p.gf_struct};

    // Initialize the result containers
    G_tau    = block_gf<imtime>{{p.beta, Fermion, p.n_tau}, p.gf_struct};
    G_iw     = G0_iw;
    Sigma_iw = G0_iw;
  }

  // -------------------------------------------------------------------------------

  void solver_core::solve(solve_params_t const &solve_params) {

    last_solve_params = solve_params;

    // http://patorjk.com/software/taag/#p=testall&f=Calvin%20S&t=TRIQS%20inchworm%0A
    if (world.rank() == 0)
      std::cout << "\n"
                   "╔╦╗╦═╗╦╔═╗ ╔═╗  ┬┌┐┌┌─┐┬ ┬┬ ┬┌─┐┬─┐┌┬┐\n"
                   " ║ ╠╦╝║║═╬╗╚═╗  │││││  ├─┤││││ │├┬┘│││\n"
                   " ╩ ╩╚═╩╚═╝╚╚═╝  ┴┘└┘└─┘┴ ┴└┴┘└─┘┴└─┴ ┴\n";

    // determine basis of operators to use
    fundamental_operator_set fops;
    for (auto const &bl : gf_struct) {
      for (auto const &a : bl.second) { fops.insert(bl.first, a); }
    }

    // setup the linear index map
    std::map<std::pair<int, int>, int> linindex;
    int block_index = 0;
    for (auto const &bl : gf_struct) {
      int inner_index = 0;
      for (auto const &a : bl.second) {
        linindex[std::make_pair(block_index, inner_index)] = fops[{bl.first, a}];
        inner_index++;
      }
      block_index++;
    }

    // Make list of block sizes
    std::vector<int> n_inner;
    for (auto const &bl : gf_struct) { n_inner.push_back(bl.second.size()); }

    // ==== Compute Delta from G0_iw ====

    auto G0_iw_inv = map([](gf_const_view<imfreq> x) { return triqs::gfs::inverse(x); }, _G0_iw);
    auto Delta_iw  = G0_iw_inv;

    for (auto &Delta_iw_bl : Delta_iw)
      for (auto const &iw : Delta_iw[0].mesh()) Delta_iw_bl[iw] = iw - Delta_iw_bl[iw];

    // Assert hermiticity of the given Weiss field
    if (!is_gf_hermitian(G0_iw)) TRIQS_RUNTIME_ERROR << "Please make sure that G0_iw fullfills the hermiticity relation G_ij[iw] = G_ji[-iw]*";

    // Merge constr_params and solve_params
    params_t params(constr_params, solve_params);

    // Reset the results
    container_set::operator=(container_set{});

    // Construct the generic Monte-Carlo solver
    triqs::mc_tools::mc_generic<mc_weight_t> mc(params.random_name, params.random_seed, params.verbosity);

    // test
    h_diag = {_h_loc, fops, params.quantum_numbers};

    // Capture random number generator
    auto &rng = mc.get_rng();

    // Create Monte-Carlo configuration
    qmc_config_t qmc_config{params,h_diag,linindex,_Delta_tau,n_inner};

    mc.add_move(moves::insert{qmc_config, rng}, "insert move");

    // Register all measurements
    if (params.measure_simple) mc.add_measure(measures::simple{params, qmc_config, result_set()}, "simple measure");

    // Perform QMC run and collect results
    mc.warmup_and_accumulate(params.n_warmup_cycles, params.n_cycles, params.length_cycle, triqs::utility::clock_callback(params.max_time));
    mc.collect_results(world);

    // Post Processing
    if (params.post_process) { post_process(params); }
  }

  // -------------------------------------------------------------------------------

  void solver_core::post_process(params_t const &p) {

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
    h5_write(grp, "G0_iw", s.G0_iw);
  }

  solver_core solver_core::h5_read_construct(triqs::h5::group h5group, std::string subgroup_name) {
    auto grp           = h5group.open_group(subgroup_name);
    auto constr_params = h5_read<constr_params_t>(grp, "constr_params");
    auto s             = solver_core{constr_params};
    h5_read(grp, "", s.result_set());
    h5_read(grp, "last_solve_params", s.last_solve_params);
    h5_read(grp, "G0_iw", s.G0_iw);
    return s;
  }

} // namespace inchworm
