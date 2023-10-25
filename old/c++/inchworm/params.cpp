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
#include "./params.hpp"

namespace inchworm {

  void h5_write(h5::group h5group, std::string subgroup_name, constr_params_t const &cp) {
    auto grp = h5group.create_group(subgroup_name);
    h5_write(grp, "n_tau", cp.n_tau);
    h5_write(grp, "n_tau_inch", cp.n_tau_inch);
    h5_write(grp, "n_tau_green", cp.n_tau_green);
    h5_write(grp, "n_iw", cp.n_iw);
    h5_write(grp, "beta", cp.beta);
    h5_write(grp, "gf_struct", cp.gf_struct);
  }

  void h5_read(h5::group h5group, std::string subgroup_name, constr_params_t &cp) {
    auto grp = h5group.open_group(subgroup_name);
    h5_read(grp, "n_tau", cp.n_tau);
    h5_read(grp, "n_tau_inch", cp.n_tau_inch);
    h5_read(grp, "n_tau_green", cp.n_tau_green);
    h5_read(grp, "n_iw", cp.n_iw);
    h5_read(grp, "beta", cp.beta);
    h5_read(grp, "gf_struct", cp.gf_struct);
  }

  void h5_write(h5::group h5group, std::string subgroup_name, solve_params_t const &sp) {
    auto grp = h5group.create_group(subgroup_name);
    h5_write(grp, "h_imp", sp.h_imp);
    h5_write(grp, "partition_method", sp.partition_method);
    h5_write(grp, "quantum_numbers", sp.quantum_numbers);
    h5_write(grp, "n_bath_sites_ED", sp.n_bath_sites_ED);
    h5_write(grp, "max_prob_zeroth_order", sp.max_prob_zeroth_order);
    h5_write(grp, "n_cycles", sp.n_cycles);
    h5_write(grp, "length_cycle", sp.length_cycle);
    h5_write(grp, "n_warmup_cycles", sp.n_warmup_cycles);
    h5_write(grp, "n_callibration_cycles", sp.n_callibration_cycles);
    h5_write(grp, "max_order", sp.max_order);
    h5_write(grp, "random_seed", sp.random_seed);
    h5_write(grp, "random_name", sp.random_name);
    h5_write(grp, "use_double_insertion", sp.use_double_insertion);
    h5_write(grp, "max_time", sp.max_time);
    h5_write(grp, "verbosity", sp.verbosity);
    h5_write(grp, "measure_average_sign", sp.measure_average_sign);
    h5_write(grp, "measure_average_order", sp.measure_average_order);
    h5_write(grp, "measure_order_histogram", sp.measure_order_histogram);
    h5_write(grp, "measure_frame_by_order", sp.measure_frame_by_order);
    h5_write(grp, "post_process", sp.post_process);
  }

  void h5_read(h5::group h5group, std::string subgroup_name, solve_params_t &sp) {
    auto grp = h5group.open_group(subgroup_name);
    // Take care! Do not read random_seed and verbosity as they should be different based on mpi rank
    h5_read(grp, "h_imp", sp.h_imp);
    h5_read(grp, "partition_method", sp.partition_method);
    h5_read(grp, "quantum_numbers", sp.quantum_numbers);
    h5::try_read(grp, "n_bath_sites_ED", sp.n_bath_sites_ED);
    h5::try_read(grp, "max_prob_zeroth_order", sp.max_prob_zeroth_order);
    h5_read(grp, "n_cycles", sp.n_cycles);
    h5_read(grp, "length_cycle", sp.length_cycle);
    h5_read(grp, "n_warmup_cycles", sp.n_warmup_cycles);
    h5::try_read(grp, "n_callibration_cycles", sp.n_callibration_cycles);
    h5::try_read(grp, "max_order", sp.max_order);
    h5_read(grp, "random_name", sp.random_name);
    h5_read(grp, "use_double_insertion", sp.use_double_insertion);
    h5_read(grp, "max_time", sp.max_time);
    h5::try_read(grp, "measure_average_sign", sp.measure_average_sign);
    h5_read(grp, "measure_average_order", sp.measure_average_order);
    h5_read(grp, "measure_order_histogram", sp.measure_order_histogram);
    h5_read(grp, "measure_frame_by_order", sp.measure_frame_by_order);
    h5_read(grp, "post_process", sp.post_process);
  }

} // namespace inchworm
