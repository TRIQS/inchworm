#include "gen_op_times.hpp"

#include "./../distributions.hpp"

namespace inchworm::moves {

  std::array<double, 2> gen_op_times(triqs::mc_tools::random_generator &rng, fop_t const &d, fop_t const &d_dag, config_t const &config,
                                     double tau_split, double tau_max) {
    EXPECTS(d.bl == d_dag.bl);

    if (tau_split == 0.0) { // ====== CTHyb sampling with bare propagator ======

      return {rng(tau_max), rng(tau_max)};

    } else { // ============ Inchworm Sampling ===============

      if (config.size() == 0) { // ------- Empty Config -------

        // Make sure that we choose tau values on seperate sides of the split points at zero and tau_split
        double dtau = tau_max - tau_split;

        if (rng(2))
          return {rng(tau_split), tau_split + rng(dtau)};
        else
          return {tau_split + rng(dtau), rng(tau_split)};

      } else { // ------ Finite size config ------

        return {get_close_time(rng, d, config.split_times, tau_max), get_close_time(rng, d_dag, config.split_times, tau_max)};
      }
    }
  }

  double get_time_prop_prob(fop_t const &d, fop_t const &d_dag, config_t const &config, double tau_split, double tau_max) {
    EXPECTS(d.bl == d_dag.bl);

    if (tau_split == 0.0) { // ----- CTHyb sampling with bare propagator

      return 1.0 / tau_max / tau_max;

    } else { // ----- Inchworm Sampling

      if (config.size() == 0) { // Account for special treatment of empty config

        // Account for insertion around tau_split and zero
        double dtau = tau_max - tau_split;

        return 1.0 / dtau / tau_split / 2.0;

      } else { // ------ Finite block size ------

        return get_prob(d, config.split_times, tau_max) * get_prob(d_dag, config.split_times, tau_max);
      }
    }
  }

} // namespace inchworm::moves
