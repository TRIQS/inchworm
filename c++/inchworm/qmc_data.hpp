#pragma once

#include "types.hpp"

namespace inchworm {

  struct config_t {
    std::vector<fop_t> d_list, d_dag_list; // list of d/d_dag not time ordered, but different
    int size() const { return d_list.size(); }
    bool try_insert(fop_t const &ddag, fop_t const &d);
    bool try_erase(int i_dag, int i);
    bool try_double_insert(fop_t const &d_dag1, fop_t const &d1, fop_t const &d_dag2, fop_t const &d2);
    bool try_double_erase(int i_dag, int i, int j_dag, int j);
  };

  struct weights_t {
    scalar_t imp; // Impurity weight - Frobenius norm of the current (propagator or green function) frame
    scalar_t hyb; // Hybridization weight - Value of the determinant in cthyb or its equivalent for the inchworm
  };

  /// The Monte-Carlo Configuration structure
  struct qmc_data_t {
    // last accepted paraemeters
    config_t config   = {};       // last accepted configuration of d and d_dag
    weights_t weights = {1., 1.}; // weight values of the last accepted configuration
    int sign          = 1;        // sign of the last accepted configuration
    frame_t frame;                // the configuration (propagator or green function) frame
  };

  // structure to calculate hybridization function for tau, tau_dag, and orbital (linear) indices.
  struct hyb_adaptor_t {
    h_tau_t const &hyb_tau;

    hyb_adaptor_t(h_tau_t const &hyb_tau) : hyb_tau(hyb_tau) {}

    // function to link
    scalar_t operator()(fop_t const &cdag, fop_t const &c) const {

      if (cdag.bl != c.bl) return 0.;
      double dtau = cdag.tau - c.tau;

      if (dtau >= 0.) {
        return hyb_tau[c.bl](dtau)(cdag.idx, c.idx);
      } else {
        return -hyb_tau[c.bl](hyb_tau[c.bl].domain().beta + dtau)(cdag.idx, c.idx);
      }
    }
  };

  // static parameters of the Monte Carlo simulation
  struct qmc_params_t {

    // same for every inch step:
    hyb_adaptor_t hyb_adaptor;
    atom_diag const &ad_imp;      // Diagonalization of the atomic problem
    bool use_bare_propagator;     // True only for the first iteration of the inchworm calculation
    MODE mode;                    // The sampling mode, either PROPAGATOR or GREENFUNCTION
    std::optional<int> max_order; // The maximum perturbation order [optional]

    // updated at every inch step:
    u_tau_t const &u_tau; // The propagator

    // different at every inch step:
    double tau_max;   // similar to beta, but configuration here does not always goes up to beta. 0 < tau_max <= beta
    double tau_split; // in the inchworm, this should be the tau_max of the previous inching. 0 < tau_split <= tau_max

    qmc_params_t(h_tau_t const &hyb_tau, atom_diag const &ad_imp, u_tau_t const &u_tau, double tau_max, double tau_split, bool use_bare_propagator,
                 MODE mode, std::optional<int> max_order)
       : hyb_adaptor(hyb_tau),
         ad_imp(ad_imp),
         use_bare_propagator(use_bare_propagator),
         mode(mode),
         max_order(max_order),
         u_tau(u_tau),
         tau_max(tau_max),
         tau_split(tau_split) {}
  };
} // namespace inchworm
