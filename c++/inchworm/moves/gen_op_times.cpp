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
        double dtau             = tau_max - tau_split;
        auto [d_tau, d_dag_tau] = rng(2) ? get_close_smaller_and_larger_time(rng, d, d_dag, tau_split, tau_split, dtau, tau_max) :
                                           get_close_smaller_and_larger_time(rng, d, d_dag, 0.0, dtau, tau_split, tau_max);

        return {d_tau, d_dag_tau};

      //} else if (config.size(d.bl) == 0) { // ------ Empty block ------

        //return get_close_smaller_and_larger_time(rng, d, d_dag, config.split_times, tau_max);

      } else { // ------ Finite block size ------

	return {get_close_time(rng, d, config.split_times, tau_max), get_close_time(rng, d_dag, config.split_times, tau_max)};

        // ------- OLD ----------
        //return {rng(tau_max), rng(tau_max)};

        // ------- NEW ----------
        // We have two options, either split-point based insertion or operator-based insertion
        //if (rng(2)) { // Operator-based insertion
                      ////FIXME Should we really weigh this with 50 percent? This path probably has significantly lower acceptance rate
	  //return {get_close_time(rng, d, config.d_dag_bl_list[d.bl], tau_max), get_close_time(rng, d_dag, config.d_bl_list[d_dag.bl], tau_max)};
        //} else { // Split-point based insertion
	  //return get_close_smaller_and_larger_time(rng, d, d_dag, config.split_times, tau_max);
        //}
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
        return 0.5
           * (get_prob_smaller_and_larger_time(d, d_dag, tau_split, tau_split, dtau, tau_max)
              + get_prob_smaller_and_larger_time(d, d_dag, 0.0, dtau, tau_split, tau_max));

      //} else if (config.size(d.bl) == 0) {
        //return get_prob_smaller_and_larger_time(d, d_dag, config.split_times, tau_max);

      } else { // ------ Finite block size ------

	return get_prob(d, config.split_times, tau_max) * get_prob(d_dag, config.split_times, tau_max);

        // ------- OLD ----------
        //return 1.0 / tau_max / tau_max;

        // ------- NEW ----------
        // We have two insertion options, either split-point based insertion or operator-based insertion
        //return 0.5
           //* (get_prob(d, config.d_dag_bl_list[d.bl], tau_max) * get_prob(d_dag, config.d_bl_list[d_dag.bl], tau_max)
              //+ get_prob_smaller_and_larger_time(d, d_dag, config.split_times, tau_max));

        //return get_prob(d, config.d_dag_bl_list[d.bl], tau_max) * get_prob(d_dag, config.d_bl_list[d_dag.bl], tau_max);
	//return get_prob_smaller_and_larger_time(d, d_dag, config.split_times, tau_max);
      }
    }
  }

}
