#pragma once

#include "../types.hpp"
#include "../solver_core.hpp"
#include "../config.hpp"

#include <triqs/mc_tools/random_generator.hpp>

namespace inchworm::moves {

  enum class KIND { Single, Double, DoubleEqBl };

  /// A simple Monte-Carlo move
  struct base_move {

    /// Attempt vertex insertion
    scalar_t attempt();

    /// Accept vertex insertion
    virtual scalar_t accept();

    /// Reject vertex insertion
    void reject() {}

    /// Constructor
    base_move(config_t &config, frame_t &frame, qmc_params_t const &params, solver_core const &solver, triqs::mc_tools::random_generator &rng)
       : config(config), frame(frame), params(params), solver(solver), rng(rng) {}

    /// Destructor
    virtual ~base_move() = default;

    /// The reweighting cutoff order
    inline static int reweighting_cutoff = 0;

    /// The reweighting coefficients
    inline static std::vector<double> reweighting_coeffs = {};

    /// Switch for enabling / disabling the tau diff statistic gathering
    inline static bool gather_tau_diff_stat = false;

    protected:
    /// The function to update the configuration
    virtual scalar_t try_config_update(config_t &) = 0;

    /// The function to update the configuration
    virtual std::string name() const = 0;

    /// The Monte-Carlo configuration
    config_t &config;

    /// The proposed Monte-Carlo configuration
    config_t prop_config = config;

    /// The propagator or Green function frame
    frame_t &frame;

    /// The proposed propagator or Green function frame
    frame_t prop_frame = frame;

    /// The Monte-Carlo parameters
    qmc_params_t const &params;

    /// The solver object
    solver_core const &solver;

    /// The random number generator
    triqs::mc_tools::random_generator &rng;

    /// The green function structure
    gf_struct_t const &gf_struct = solver.constr_params.gf_struct;
  };

} // namespace inchworm::moves
