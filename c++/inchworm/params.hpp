/*******************************************************************************
 *
 * inchworm: A TRIQS based impurity solver
 *
 * Copyright (c) 2019 The Simons foundation
 *   authors: Maxime Charlebois
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
#pragma once

#include "types.hpp"

namespace inchworm {

  /// The parameters for the solver construction
  struct constr_params_t {

    /// Number of tau points for the hybridization function
    int n_tau = 101;

    /// Number of tau points for the propagator
    int n_tau_inch = 11;

    /// Number of tau points for the Green function
    int n_tau_green = 101;

    /// Number of Matsubara frequencies
    int n_iw = 5;

    /// Inverse temperature
    double beta;

    /// Block structure of the gf
    gf_struct_t gf_struct;

    /// Number of block indeces for the Green function
    int n_blocks() const { return gf_struct.size(); }

    /// Names of block indeces for the Green function
    auto block_names() const {
      std::vector<std::string> v;
      for (auto const &bl : gf_struct) v.push_back(bl.first);
      return v;
    }

    /// Write constr_params_t to hdf5
    friend void h5_write(h5::group h5group, std::string subgroup_name, constr_params_t const &cp);

    /// Read constr_params_t from hdf5
    friend void h5_read(h5::group h5group, std::string subgroup_name, constr_params_t &cp);
  };

  /// The parameters for the solve function
  struct solve_params_t {

    // ----------- System Specific -----------

    /// Impurity Hamiltonian
    many_body_operator h_imp;

    /// Partition method
    /// type: str
    std::string partition_method = "automatic";

    /// Quantum numbers
    /// type: list(Operator)
    /// default: [Total Particle Number]
    std::vector<many_body_op_t> quantum_numbers = {};

    // ----------- QMC Specific -----------

    /// Number of MC cycles
    int n_cycles;

    /// Length of a MC cycles
    int length_cycle = 50;

    /// Number of warmup cycles
    int n_warmup_cycles = 5000;

    /// Random seed of the random generator
    int random_seed = 34789 + 928374 * mpi::communicator().rank();

    /// Name of the random generator
    std::string random_name = "";

    /// Use double insertion
    int use_double_insertion = true;

    /// Maximum running time in seconds (-1 : no limit)
    int max_time = -1;

    /// Verbosity
    int verbosity = mpi::communicator().rank() == 0 ? 1 : 0;

    // ----------- Measurements -----------

    /// Measure the average perturbation order
    bool measure_average_order = true;

    /// Measure the average perturbation order
    bool measure_order_histogram = false;

    /// Measure the average perturbation order
    bool measure_frame_by_order = false;

    /// Perform post processing
    bool post_process = true;

    /// Write constr_params_t to hdf5
    friend void h5_write(h5::group h5group, std::string subgroup_name, solve_params_t const &sp);

    /// Read constr_params_t from hdf5
    friend void h5_read(h5::group h5group, std::string subgroup_name, solve_params_t &sp);
  };

  /// A struct combining both constr_params_t and solve_params_t
  struct params_t : constr_params_t, solve_params_t {
    params_t(constr_params_t const &constr_params_, solve_params_t const &solve_params_)
       : constr_params_t(constr_params_), solve_params_t(solve_params_) {}
  };

} // namespace inchworm
