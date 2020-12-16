#pragma once

#include "../types.hpp"
#include "../qmc_data.hpp"

#include <triqs/mc_tools/random_generator.hpp>

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
    base_move(qmc_data_t &data, gf_struct_t const &gf_struct, qmc_params_t const &qmc_params, triqs::mc_tools::random_generator &rng);

    /// Destructor
    virtual ~base_move() = default;

    /// The cutoff for the weighting function
    inline static long reweighting_cutoff = 0;

    protected:
    /// The function to update the configuration
    virtual scalar_t try_config_update(config_t &) = 0;

    /// The function to update the configuration
    virtual std::string name() const = 0;

    /// The vector of all creation operators
    std::vector<fop_t> all_d_ops;

    /// The vector of all annihilation operators
    std::vector<fop_t> all_d_dag_ops;

    /// The Monte-Carlo configuration
    qmc_data_t &data;

    /// The Monte-Carlo configuration
    qmc_data_t prop_data;

    /// The Monte-Carlo parameters
    qmc_params_t const &params;

    /// The random number generator
    triqs::mc_tools::random_generator &rng;

    /// the green function structure
    gf_struct_t const &gf_struct;
  };

} // namespace inchworm::moves
