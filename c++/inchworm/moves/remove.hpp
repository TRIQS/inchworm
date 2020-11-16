#pragma once
#include "./base_move.hpp"

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
    scalar_t try_move(config_t & config) override;
  };

  struct double_remove : public base_move {
    using base_move::base_move;
    inline std::string name() const override { return "Double Remove"; };

    /**
     * Try to remove two (cdag, c) operator pairs from the configuration
     *
     * @param config The configuration to update
     * @return On success, return the ratio between the inverse move probability
     *         and the move probability.
     *         On failure, return 0.0
     */
    scalar_t try_move(config_t & config) override;
  };

} // namespace inchworm::moves
