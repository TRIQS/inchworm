#include "./insert.hpp"
#include "./../torus.hpp"
#include "./../distributions.hpp"

namespace inchworm::moves {

  scalar_t insert::try_config_update(config_t &config) {
    long n_bl = gf_struct.size();

    long bl      = rng(n_bl);
    long bl_size = gf_struct[bl].second;

    auto d     = solver.all_d_ops[bl][rng(bl_size)];
    auto d_dag = solver.all_d_dag_ops[bl][rng(bl_size)];

    scalar_t t_ratio;
    if (params.tau_split == 0.0) { // ----- CTHyb sampling with bare propagator

      // Choose random times uniformly
      d.tau     = rng(params.tau_max);
      d_dag.tau = rng(params.tau_max);

      long new_nop_bl = config.size(bl) + 1;
      t_ratio         = std::pow(params.tau_max * bl_size / new_nop_bl, 2);

    } else { // ----- Inchworm Sampling

      if (config.size() == 0) {
        // Special case of empty configuration for inchworm
        // Make sure that we choose tau values on seperate
        // sides of the split point

        auto dtau = params.tau_max - params.tau_split;
        d.tau     = rng(params.tau_split);
        d_dag.tau = params.tau_split + rng(dtau);
        if (rng(2)) std::swap(d.tau, d_dag.tau);
        t_ratio = 2.0 * params.tau_split * dtau * bl_size * bl_size;
      } else {

        // Choose the d and d_dag to insert next to
        // FIXME Make this Block Specific!!!
        // Requires special treatment of (config.size(bl) == 0) however ..

        //long old_nop_bl = config.size(bl);
        //double d_tau_ref     = config.d_bl_list[bl][rng(old_nop_bl)].tau;
        //double d_dag_tau_ref = config.d_dag_bl_list[bl][rng(old_nop_bl)].tau;
        //t_ratio = std::pow(double(bl_size) * old_nop_bl / new_nop_bl, 2) / prob_d_tau / prob_d_dag_tau;

        auto [d_tau, prob_d_tau]         = get_close_time_and_prob(rng, config.d_dag_list, params.tau_max);
        auto [d_dag_tau, prob_d_dag_tau] = get_close_time_and_prob(rng, config.d_list, params.tau_max);

        d.tau     = d_tau;
        d_dag.tau = d_dag_tau;

        long new_nop_bl = config.size(bl) + 1;
        t_ratio         = std::pow(double(bl_size) / new_nop_bl, 2) / prob_d_tau / prob_d_dag_tau;
      }
    }

    if (not config.try_insert(d_dag, d)) return 0.0;

    return t_ratio;
  }

  scalar_t double_insert::try_config_update(config_t &config) {
    long n_bl = gf_struct.size();

    long bl1, bl2;
    if (equal_blocks) {
      bl1 = rng(n_bl);
      bl2 = bl1;
    } else {
      bl1 = rng(n_bl);
      bl2 = rng(n_bl);
    }

    long bl1_size = gf_struct[bl1].second;
    long bl2_size = gf_struct[bl2].second;

    auto d1     = solver.all_d_ops[bl1][rng(bl1_size)];
    auto d1_dag = solver.all_d_dag_ops[bl1][rng(bl1_size)];
    auto d2     = solver.all_d_ops[bl2][rng(bl2_size)];
    auto d2_dag = solver.all_d_dag_ops[bl2][rng(bl2_size)];

    d1.tau     = rng(params.tau_max);
    d1_dag.tau = rng(params.tau_max);
    d2.tau     = rng(params.tau_max);
    d2_dag.tau = rng(params.tau_max);

    if (not config.try_double_insert(d1_dag, d1, d2_dag, d2)) return 0;

    long new_nop_bl1 = config.d_bl_list[bl1].size();
    long new_nop_bl2 = config.d_bl_list[bl2].size();
    return std::pow(params.tau_max * params.tau_max * bl1_size * bl2_size / (new_nop_bl1 * new_nop_bl2), 2);
  }

} // namespace inchworm::moves
