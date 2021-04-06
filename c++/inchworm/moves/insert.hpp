#pragma once

#include "base_move.hpp"

#include <numeric>

namespace inchworm::moves {

  struct insert : public base_move {

    using acc_t = std::vector<nda::array<accumulator<double>, 3>>;

    insert(config_t &config, frame_t &frame, qmc_params_t const &params, solver_core const &solver, triqs::mc_tools::random_generator &rng,
           acc_t &tau_diff_stat)
       : base_move(config, frame, params, solver, rng), tau_diff_stat(tau_diff_stat) {}

    inline std::string name() const override { return "Insert"; };

    /**
     * Try to insert a single (cdag, c) operator pair into the configuration
     *
     * @param config The configuration to update
     * @return On success, return the ratio between the inverse move probability
     *         and the move probability.
     *         On failure, return 0.0
     */
    scalar_t try_config_update(config_t &config) override;

    scalar_t accept() override;

    double beta = solver.constr_params.beta;

    fop_t last_d     = {};
    fop_t last_d_dag = {};

    inline static bool gather_tau_diff_stat = false;
    acc_t &tau_diff_stat;
  };

  struct double_insert : public base_move {

    bool equal_blocks;

    double_insert(config_t &config, frame_t &frame, qmc_params_t const &params, solver_core const &solver, triqs::mc_tools::random_generator &rng,
                  bool equal_blocks)
       : base_move(config, frame, params, solver, rng), equal_blocks(equal_blocks) {}

    inline std::string name() const override {
      if (equal_blocks) return "Double Insert (equal_blocks)";
      return "Double Insert";
    };

    /**
     * Try to insert two (cdag, c) operator pairs into the configuration
     *
     * @param config The configuration to update
     * @return On success, return the ratio between the inverse move probability
     *         and the move probability.
     *         On failure, return 0.0
     */
    scalar_t try_config_update(config_t &config) override;
  };

} // namespace inchworm::moves
