#pragma once

#include "base_move.hpp"

namespace inchworm::moves {

  struct insert : public base_move {
    using base_move::base_move;
    inline std::string name() const override { return "Insert"; };

    /**
     * Try to insert a single (cdag, c) operator pair into the configuration
     *
     * @param config The configuration to update
     * @return On success, return the ratio between the inverse move probability
     *         and the move probability.
     *         On failure, return 0.0
     */
    scalar_t try_config_update(config_t & config) override;
  };

  struct double_insert : public base_move {

    bool equal_blocks;

    double_insert(qmc_data_t &data, gf_struct_t const &gf_struct, qmc_params_t const &qmc_params, triqs::mc_tools::random_generator &rng,
                  bool equal_blocks)
       : base_move(data, gf_struct, qmc_params, rng), equal_blocks(equal_blocks) {}

    inline std::string name() const override { return "Double Insert"; };

    /**
     * Try to insert two (cdag, c) operator pairs into the configuration
     *
     * @param config The configuration to update
     * @return On success, return the ratio between the inverse move probability
     *         and the move probability.
     *         On failure, return 0.0
     */
    scalar_t try_config_update(config_t & config) override;
  };

} // namespace inchworm::moves
