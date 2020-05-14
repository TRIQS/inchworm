#pragma once
#include "./container_set.hpp"
#include "./params.hpp"
#include "./types.hpp"
#include "./mc/impurity_product.hpp"
#include "./u_frame.hpp"
#include "./util.hpp"

namespace inchworm {

  /// The Solver class
  class solver_core : public container_set {

    private:
    double beta;           // inverse temperature
    atom_diag h_diag;      // diagonalization of the local problem
    gf_struct_t gf_struct; // Block structure of the Green function FIXME
    many_body_op_t _h_loc; // The local Hamiltonian = h_int + h0
    std::map<int, std::pair<int, int>> map_lin_idx_to_block_inner;
    fundamental_operator_set fops;
    u_frame_t u_frame_bare;

    //mpi::communicator _comm;   // define the communicator, here MPI_COMM_WORLD
    int _solve_status; // Status of the solve upon exit: 0 for clean termination, > 0 otherwise.

    // Single-particle Green's function containers
    // std::vector<matrix<dcomplex>> Delta_infty_vec; // Quadratic instantaneous part of G0_iw

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

    void init(solve_params_t const &solve_params);

    single_step_results_t single_step(solve_params_t const &solve_params, double tau_split, double tau_max, bool use_bare_propagator);
    single_step_results_t solve_single_step(solve_params_t const &solve_params, bool use_bare_propagator, double tau_split, double tau_max);
    single_step_results_t solve_self_consistently(solve_params_t const &solve_params, u_tau_t const &u_tau_, double tau_split, double tau_max);

    void solve_inchworm(solve_params_t const &solve_params);

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

    static std::string hdf5_scheme() { return "INCHWORM_SolverCore"; }

    // Function that writes a solver object to hdf5 file
    friend void h5_write(triqs::h5::group h5group, std::string subgroup_name, solver_core const &s);

    // Function that constructs a solver object from an hdf5 file
    CPP2PY_IGNORE
    static solver_core h5_read_construct(triqs::h5::group h5group, std::string subgroup_name);
  };
} // namespace inchworm
