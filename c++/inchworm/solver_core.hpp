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
#pragma once
#include "./container_set.hpp"
#include "./params.hpp"
#include "./types.hpp"

namespace inchworm {

  /// The Solver class
  class solver_core : public container_set {

    private:
    double beta;           // inverse temperature
    atom_diag h_diag;      // diagonalization of the local problem
    gf_struct_t gf_struct; // Block structure of the Green function FIXME
    many_body_op_t _h_loc; // The local Hamiltonian = h_int + h0

    //mpi::communicator _comm;   // define the communicator, here MPI_COMM_WORLD
    int _solve_status;         // Status of the solve upon exit: 0 for clean termination, > 0 otherwise.

    // Single-particle Green's function containers
    g_iw_t _G0_iw;                                 // Non-interacting Matsubara Green's function
    g_tau_t _Delta_tau;                            // Imaginary-time Hybridization function
    std::vector<matrix<dcomplex>> Delta_infty_vec; // Quadratic instantaneous part of G0_iw

    // Mpi Communicator
    mpi::communicator world;

    // Return reference to container_set
    container_set &result_set() { return static_cast<container_set &>(*this); }
    container_set const &result_set() const { return static_cast<container_set const &>(*this); }

    // Function to perform the post-processing steps
    void post_process(params_t const &p);

    public:
    /**
     * Construct a INCHWORM solver
     *
     * @param construct_parameters Set of parameters specific to the INCHWORM solver
     */
    CPP2PY_ARG_AS_DICT
    solver_core(constr_params_t const &constr_params_);

    // Delete assignement operator because of const members
    solver_core(solver_core const &p) = default;
    solver_core(solver_core &&p)      = default;
    solver_core &operator=(solver_core const &p) = delete;
    solver_core &operator=(solver_core &&p) = default;

    /**
     * Solve method that performs INCHWORM calculation
     *
     * @param solve_params_t Set of parameters specific to the INCHWORM run
     */
    CPP2PY_ARG_AS_DICT
    void solve(solve_params_t const &solve_params);

    // Struct containing the parameters relevant for the solver construction
    constr_params_t constr_params;

    // Struct containing the parameters relevant for the solve process
    std::optional<solve_params_t> last_solve_params;

    /// Noninteracting Green Function in Matsubara frequencies
    g_iw_t G0_iw;

    // Allow the user to retrigger post-processing with the last set of parameters
    void post_process() {
      if (not last_solve_params) TRIQS_RUNTIME_ERROR << "You need to run the solver once before you post-process";
      post_process({constr_params, last_solve_params.value()});
    }

    static std::string hdf5_scheme() { return "INCHWORM_SolverCore"; }

    // Function that writes a solver object to hdf5 file
    friend void h5_write(triqs::h5::group h5group, std::string subgroup_name, solver_core const &s);

    // Function that constructs a solver object from an hdf5 file
    CPP2PY_IGNORE
    static solver_core h5_read_construct(triqs::h5::group h5group, std::string subgroup_name);
  };
} // namespace inchworm
