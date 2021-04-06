#pragma once

#include "base_move.hpp"

namespace inchworm::moves {

  struct remove : public base_move {
    using base_move::base_move;
    inline std::string name() const override { return "Remove"; };

    /**
     * Try to remove a single (cdag, c) operator pair from the configuration
     *
     * @param config The configuration to update
     * @return On success, return the ratio between the inverse move probability
     *         and the move probability.
     *         On failure, return 0.0
     */
    scalar_t try_config_update(config_t & config) override;
  };

  struct double_remove : public base_move {

    bool equal_blocks;

    double_remove(qmc_data_t &data, qmc_params_t const &params, solver_core const &solver, triqs::mc_tools::random_generator &rng, bool equal_blocks)
       : base_move(data, params, solver, rng), equal_blocks(equal_blocks) {}

    inline std::string name() const override {
      if (equal_blocks) return "Double Remove (equal_blocks)";
      return "Double Remove";
    };

    /**
     * Try to remove two (cdag, c) operator pairs from the configuration
     *
     * @param config The configuration to update
     * @return On success, return the ratio between the inverse move probability
     *         and the move probability.
     *         On failure, return 0.0
     */
    scalar_t try_config_update(config_t & config) override;
  };

} // namespace inchworm::moves
