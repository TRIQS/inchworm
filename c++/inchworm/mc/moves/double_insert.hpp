#pragma once
#include "../qmc_config_data.hpp"
#include "../../diagram/print.hpp"
#include "../../u_frame.hpp"

#include <triqs/mc_tools/random_generator.hpp>
#include <inchworm/diagram/inclusion_exclusion.hpp>
#include <inchworm/diagram/proper_enum.hpp>

namespace inchworm::moves {

  /// A simple Monte-Carlo move
  struct double_insert {

    /// Attempt vertex double insertion
    scalar_t attempt();

    /// Accept vertex double insertion
    scalar_t accept();

    /// Reject vertex double insertion
    void reject() {}

    /// Constructor
    double_insert(qmc_config_data_t &data, params_t const &params, qmc_params_t const &qmc_params, triqs::mc_tools::random_generator &rng)
       : data(data), params(qmc_params), rng(rng), gf_struct(params.gf_struct) {
      prop_g_frame = make_frame(params.gf_struct);

      for (auto const &op : qmc_params.ad_imp.get_fops()) {
        auto bl_name = std::get<std::string>(op.index[0]);
        auto idx     = std::get<long>(op.index[1]);

        // Determine the number of the bl_name in gf_struct
        auto it = std::find_if(gf_struct.cbegin(), gf_struct.cend(), [&](auto &&x) { return x.first == bl_name; });
        long bl = std::distance(gf_struct.cbegin(), it);

        all_d_ops.push_back({0.0, false, op.linear_index, bl, idx});
        all_d_dag_ops.push_back({0.0, true, op.linear_index, bl, idx});
      }
    }

    private:
    /// The vector of all valid creation operators
    std::vector<fop_t> all_d_ops;

    /// The vector of all annihilation operators
    std::vector<fop_t> all_d_dag_ops;

    /// The Monte-Carlo configuration
    qmc_config_data_t &data;

    /// The Monte-Carlo parameters
    qmc_params_t const &params;

    /// The random number generator
    triqs::mc_tools::random_generator &rng;

    /// d/d_dag lists of proposed double_insert
    config_t prop_config; // proposed configuration of d and d_dag

    /// weights of the proposed configuration
    weights_t prop_weights;

    /// container of the calculated time frame of the propagator
    u_partial_t prop_u_partial;

    /// the green function structure
    gf_struct_t const &gf_struct;

    /// green function container to accumulate into
    frame_t prop_g_frame;

    int prop_sign;
  };

} // namespace inchworm::moves
