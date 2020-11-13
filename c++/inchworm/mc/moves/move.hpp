#pragma once
#include "../qmc_config_data.hpp"
#include "../../diagram/print.hpp"
#include "../../u_frame.hpp"

#include <triqs/mc_tools/random_generator.hpp>
#include <inchworm/diagram/inclusion_exclusion.hpp>
#include <inchworm/diagram/proper_enum.hpp>

namespace inchworm::moves {

  /// A simple Monte-Carlo move
  struct base_move {

    /// Attempt vertex insertion
    scalar_t attempt();

    /// Accept vertex insertion
    scalar_t accept();

    /// Reject vertex insertion
    void reject() {}

    /// Constructor FIXME params naming
    base_move(qmc_config_data_t &data, gf_struct_t const &gf_struct, qmc_params_t const &qmc_params, triqs::mc_tools::random_generator &rng);

    /// Destructor
    virtual ~base_move() = default;

    protected:
    /// The function to update the configuration
    virtual scalar_t try_move(config_t &) = 0;

    /// The function to update the configuration
    virtual std::string name() const = 0;

    /// The vector of all creation operators
    std::vector<fop_t> all_d_ops;

    /// The vector of all annihilation operators
    std::vector<fop_t> all_d_dag_ops;

    /// The Monte-Carlo configuration
    qmc_config_data_t &data;

    /// The Monte-Carlo configuration
    qmc_config_data_t prop_data;

    /// The Monte-Carlo parameters
    qmc_params_t const &params;

    /// The random number generator
    triqs::mc_tools::random_generator &rng;

    /// the green function structure
    gf_struct_t const &gf_struct;
  };

  struct insert : public base_move {
    using base_move::base_move;
    inline std::string name() const override { return "Insert"; };
    scalar_t try_move(config_t &) override;
  };

  struct double_insert : public base_move {
    using base_move::base_move;
    inline std::string name() const override { return "Double Insert"; };
    scalar_t try_move(config_t &) override;
  };

  struct remove : public base_move {
    using base_move::base_move;
    inline std::string name() const override { return "Remove"; };
    scalar_t try_move(config_t &) override;
  };

  struct double_remove : public base_move {
    using base_move::base_move;
    inline std::string name() const override { return "Double Remove"; };
    scalar_t try_move(config_t &) override;
  };

} // namespace inchworm::moves
