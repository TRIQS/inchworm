#pragma once
#include "../u_frame.hpp"
#include "../types.hpp"
#include "../params.hpp"
#include "../diagram/diagram.hpp"
#include "./impurity_product.hpp"
#include <triqs/atom_diag/atom_diag.hpp>

namespace inchworm {

  struct weights_t {
    scalar_t loc; // atomic weight (Frobenius norm of the current propagator frame)
    scalar_t hyb; // value of the determinant in cthyb or its equivalent for the inchworm
  };

  struct config_t {
    std::vector<diagram::time_and_index_t> d_list, d_dag_list; // list of d/d_dag not time ordered, but different
    int size() const { return d_list.size(); }
    bool try_insert(double tau, int linear_index, double tau_dag, int linear_index_dag);
    bool try_erase(int i, int i_dag);
  };

  /// The Monte-Carlo Configuration structure
  struct qmc_config_data_t {

    // last accepted paraemeters
    config_t config;       // last accepted configuration of d and d_dag
    weights_t w;           // weight values of the last accepted configuration
    u_partial_t u_partial; // frame of the last accepted configuraiton: just one time frame of a propagator
    frame_t g_frame;       // green function for a specific time tau
    int sign;              // sign of the last accepted configuration
    //scalar_t normalization_cte; //

    qmc_config_data_t(gf_struct_t const & gf_struct) : w{1., 1.}, g_frame{make_frame(gf_struct)}, sign{1} {}
  };

  // structure to calculate hybridization function for tau, tau_dag, and orbital (linear) indices.
  struct hyb_adaptor_t {
    h_tau_t const hyb_tau;
    std::map<int, std::pair<int, int>> const &map_lin_idx_to_block_inner;

    hyb_adaptor_t(h_tau_t const &hyb_tau, std::map<int, std::pair<int, int>> const &map_lin_idx_to_block_inner)
       : hyb_tau(std::move(hyb_tau)), map_lin_idx_to_block_inner(std::move(map_lin_idx_to_block_inner)) {}

    // function to link
    scalar_t operator()(double tau, int li, double tau_dag, int li_dag) const {
      auto [bl, in]         = map_lin_idx_to_block_inner.at(li);
      auto [bl_dag, in_dag] = map_lin_idx_to_block_inner.at(li_dag);

      if (bl != bl_dag) return 0.; // important: there should be no finite terms of the hybridization between different [bl]ock.
      double dtau = tau_dag - tau;
      //std::printf("salut % 4.8f  % 4.8f  % 4.8f \n",dtau, tau_dag, tau);
      if (dtau >= 0.) {
        //std::printf("test1: %d %d   % 4.8e\n", (hyb_tau[bl])(dtau)(in_dag, in), in_dag, in);
        return (hyb_tau[bl])(dtau)(in_dag, in);
      } else {
        //std::printf("test2: %d %d   % 4.8e\n", -(hyb_tau[bl])(hyb_tau[bl].domain().beta + dtau)(in_dag, in), in_dag, in);
        return -(hyb_tau[bl])(hyb_tau[bl].domain().beta + dtau)(in_dag, in);
      }
    }
  };

  // static parameters of the Monte Carlo simulation
  struct qmc_params_t {

    // same for every inch step:
    hyb_adaptor_t hyb_adaptor;
    std::map<int, std::pair<int, int>> map_lin_idx_to_block_inner;
    atom_diag const &ad_imp; // Diagonalization of the atomic problem

    // updated at every inch step:
    u_tau_t const &u_tau; // THE propagator

    // different at every inch step:
    double tau_max;           // similar to beta, but configuration here does not always goes up to beta. 0 < tau_max <= beta
    double tau_split;         // in the inchworm, this should be the tau_max of the previous inching. 0 < tau_split <= tau_max
    bool use_bare_propagator; // true only for the first iteration of the inchworm calculation

    int mode;

    qmc_params_t(h_tau_t const &hyb_tau, std::map<int, std::pair<int, int>> const &map_lin_idx_to_block_inner, atom_diag const &ad_imp, u_tau_t const &u_tau,
                 double tau_max, double tau_split, bool use_bare_propagator, int mode)
       : hyb_adaptor(hyb_tau, map_lin_idx_to_block_inner),
         map_lin_idx_to_block_inner(map_lin_idx_to_block_inner),
         ad_imp(ad_imp),
         u_tau(u_tau),
         tau_max(tau_max),
         tau_split(tau_split),
         use_bare_propagator(use_bare_propagator),
         mode(mode) {}
  };
} // namespace inchworm
