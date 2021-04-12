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

      if (config.size() == 0) { // Account for special treatment of empty config

	// Account for insertion around tau_split and zero
	double dtau = params.tau_max - params.tau_split;
	double inverse_prop_prob = 0.5
	   * (get_prob_smaller_and_larger_time(d, d_dag, params.tau_split, params.tau_split, dtau, params.tau_max)
	      + get_prob_smaller_and_larger_time(d, d_dag, 0.0, dtau, params.tau_split, params.tau_max));
	return inverse_prop_prob / bl_size / bl_size;

      } else if (config.size(bl) == 0) {
	double inverse_prop_prob = get_prob_smaller_and_larger_time(d, d_dag, config.split_times, params.tau_max);
	return inverse_prop_prob / bl_size / bl_size;

      } else {

        // We have two insertion options, either split-point based insertion or operator-based insertion
        double inverse_prop_prob = 0.5
           * (get_prob(d, config.d_dag_bl_list[bl], params.tau_max) * get_prob(d_dag, config.d_bl_list[bl], params.tau_max)
              + get_prob_smaller_and_larger_time(d, d_dag, config.split_times, params.tau_max));
        return std::pow(double(old_nop_bl) / bl_size, 2) * inverse_prop_prob;
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
