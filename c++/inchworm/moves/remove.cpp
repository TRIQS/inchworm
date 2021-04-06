#include "./remove.hpp"
#include "./../torus.hpp"
#include "./../distributions.hpp"

namespace inchworm::moves {

  scalar_t remove::try_config_update(config_t &config) {
    long n_bl = gf_struct.size();

    long bl         = rng(n_bl);
    long old_nop_bl = config.size(bl);

    if (old_nop_bl == 0) return 0;

    long idx     = rng(old_nop_bl);
    long idx_dag = rng(old_nop_bl);

    fop_t d     = config.d_bl_list[bl][idx];
    fop_t d_dag = config.d_dag_bl_list[bl][idx];

    if (not config.try_erase(bl, idx_dag, idx)) return 0;

    long bl_size = gf_struct[bl].second;

    if (params.tau_split == 0.0) { // ----- CTHyb sampling with bare propagator

      return std::pow(old_nop_bl / (params.tau_max * bl_size), 2);

    } else { // ----- Inchworm Sampling

      // Account for special treatment of insertion into empty configuration
      if (config.size() == 0) {
        auto dtau = params.tau_max - params.tau_split;
        return 1.0 / (2.0 * params.tau_split * dtau * bl_size * bl_size);
      } else {

        double prob_d_tau     = get_prob(d, config.d_dag_list, params.tau_max);
        double prob_d_dag_tau = get_prob(d_dag, config.d_list, params.tau_max);
        return std::pow(double(old_nop_bl) / bl_size, 2) * prob_d_tau * prob_d_dag_tau;
      }
    }
  }

  scalar_t double_remove::try_config_update(config_t &config) {
    long n_bl = gf_struct.size();

    long bl1, bl2;
    if (equal_blocks) {
      bl1 = rng(n_bl);
      bl2 = bl1;
    } else {
      bl1 = rng(n_bl);
      bl2 = rng(n_bl);
    }

    long old_nop_bl1 = config.size(bl1);
    long old_nop_bl2 = config.size(bl2);

    if (old_nop_bl1 == 0 or old_nop_bl2 == 0) return 0;

    long idx1     = rng(old_nop_bl1);
    long idx1_dag = rng(old_nop_bl1);
    long idx2     = rng(old_nop_bl2);
    long idx2_dag = rng(old_nop_bl2);

    if (not config.try_double_erase(bl1, idx1_dag, idx1, bl2, idx2_dag, idx2)) return 0;

    long bl1_size = gf_struct[bl1].second;
    long bl2_size = gf_struct[bl2].second;
    return std::pow(old_nop_bl1 * old_nop_bl2 / (params.tau_max * params.tau_max * bl1_size * bl2_size), 2);
  }

} // namespace inchworm::moves
