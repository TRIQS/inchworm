#pragma once

#include "container_set.hpp"
#include "params.hpp"
#include "types.hpp"

namespace inchworm {

  /// The Solver class
  class solver_core : public container_set {

    public:
    /// The propagator in imaginary time
    u_tau_t u_tau;

    /// Diagonalization of the local problem
    atom_diag ad_imp;

    private:
    // The fundamental operator set associated with constr_params.gf_struct
    fundamental_operator_set fops;

    // Mapping of linear operator index to (block, orbital)
    std::map<int, std::pair<int, int>> map_lin_idx_to_block_inner;

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

    // solve cthyb (no split point + bare propagator):
    CPP2PY_IGNORE
    qmc_step_results_t solve_cthyb(solve_params_t const &solve_params, double tau_max);

    // self consistent solution (one step, with precalculated U(beta) from ED)
    CPP2PY_IGNORE
    qmc_step_results_t solve_self_consistently(solve_params_t const &solve_params, u_tau_t const &u_tau_, double tau_split, double tau_max);

    // Run inchworm to calculate S.u_tau
    CPP2PY_ARG_AS_DICT
    void solve_inchworm(solve_params_t const &solve_params);

    // Sample the Green function S.G_tau
    CPP2PY_ARG_AS_DICT
    void solve_green(solve_params_t const &solve_params);

    private:
    // Initialize the solver object
    void init(solve_params_t const &solve_params);

    // one Monte Carlo step calculation (common to the 3 solve scheme above):
    qmc_step_results_t qmc_step(solve_params_t const &solve_params, double tau_split, double tau_max, bool use_bare_propagator, MODE mode);

    public:
    // Struct containing the parameters relevant for the solve process
    std::optional<solve_params_t> last_solve_params;

    // Imaginary-time Hybridization function
    h_tau_t Delta_tau;

    g_iw_t G0_iw; // Non-interacting Matsubara Green's function

    // Allow the user to retrigger post-processing with the last set of parameters
    void post_process() {
      if (not last_solve_params) TRIQS_RUNTIME_ERROR << "You need to run the solver once before you post-process";
      post_process({constr_params, last_solve_params.value()});
    }

    static std::string hdf5_format() { return "INCHWORM_SolverCore"; }

    // Function that writes a solver object to hdf5 file
    friend void h5_write(h5::group h5group, std::string subgroup_name, solver_core const &s);

    // Function that constructs a solver object from an hdf5 file
    CPP2PY_IGNORE
    static solver_core h5_read_construct(h5::group h5group, std::string subgroup_name);
  };

} // namespace inchworm
