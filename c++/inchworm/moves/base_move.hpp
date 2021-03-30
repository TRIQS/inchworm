#pragma once

#include "../types.hpp"
#include "../solver_core.hpp"
#include "../config.hpp"

#include <triqs/mc_tools/random_generator.hpp>

#include <fmt/core.h>

namespace inchworm::moves {

  enum class KIND { Single, Double, DoubleEqBl };

  /// A simple Monte-Carlo move
  struct base_move {

    /// Attempt vertex insertion
    scalar_t attempt();

    /// Accept vertex insertion
    virtual scalar_t accept();

    /// Reject vertex insertion
    void reject() {
      ++reject_stats[name()][config.size()][1];
      if (print_condition()) fmt::print("Reject with reason: {}\n", rejection_reason);
      ++reject_count;
    }

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

    /// Are we accumulating yet?
    inline static bool accumulating = false;

    /// Max Moves since last empty configuration
    inline static long long max_moves_since_empty = 0;

    /// Max Expensive Moves since last empty configuration
    inline static long long max_expensive_moves_since_empty = 0;

    /// Max Accepts since last empty configuration
    inline static long long max_accepts_since_empty = 0;

    /// Rejection statistic vector
    inline static std::map<std::string, std::vector<sso_vector<long long>>> reject_stats = {};
    inline static const std::vector<std::string> stat_names = {"NMoves",    "AccRate",   "try_move",  "max_order", "quick hyb",
                                                               "quick imp", "ratio"};

    static void print_reject_stats();

    protected:
    /// The function to update the configuration
    virtual double try_config_update(config_t &) = 0;

    /// The function to update the configuration
    virtual std::string name() const = 0;

    /// The function to update the configuration
    bool print_condition() const;

    /// Print the move statistics
    void print_move_stats() const;

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

    /// Rejection Count
    inline static long long reject_count = 0;

    /// Moves since last empty configuration
    inline static long long moves_since_empty = 0;

    /// Expensive Moves since last empty configuration
    inline static long long expensive_moves_since_empty = 0;

    /// Accepts since last empty configuration
    inline static long long accepts_since_empty = 0;

    /// Reason for last rejection
    std::string rejection_reason = "";
  };

} // namespace inchworm::moves
