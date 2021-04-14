#include "./remove.hpp"
#include "./gen_op_times.hpp"

namespace inchworm::moves {

  scalar_t remove::try_config_update(config_t &config) {
    long n_bl = gf_struct.size();

    long bl         = rng(n_bl);

    if (config.size(bl) == 0) return 0;

    long idx     = rng(config.size(bl));
    long idx_dag = rng(config.size(bl));

    fop_t d     = config.d_bl_list[bl][idx];
    fop_t d_dag = config.d_dag_bl_list[bl][idx];

    double prop_prob = 1.0 / config.size(bl) / config.size(bl) / n_bl;

    if (not config.try_erase(bl, idx_dag, idx)) return 0;

    long bl_size = gf_struct[bl].second;
    double inv_prop_prob = get_time_prop_prob(d, d_dag, config, params.tau_split, params.tau_max) / bl_size / bl_size / n_bl;

    return inv_prop_prob / prop_prob;
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
