#pragma once
#include "../types.hpp"
#include "../params.hpp"
#include "../util.hpp"
#include "../diagram/diagram.hpp"
#include "./impurity_product.hpp"
#include <triqs/atom_diag/atom_diag.hpp>
//#include <triqs/atom_diag/functions.hpp>
//#include <triqs/utility/serialization.hpp>
//#include <triqs/det_manip.hpp>

namespace inchworm {

  struct weights_t {
    scalar_t loc; // atomic weight (Frobenius norm of the current propagator frame)
    scalar_t hyb; // value of the determinant in cthyb or its equivalent for the inchworm
  };

  struct config_t {
    std::vector<diagram::time_and_index_t> c_list, cdag_list; // list of c/cdag not time ordered, but different
    int size() { return c_list.size(); }
    bool try_insert(double tau, int linear_index, double tau_dag, int linear_index_dag);
    bool try_erase(int i, int i_dag);
  };

  /// The Monte-Carlo Configuration Class
  struct qmc_config_data_t {

    // last accepted paraemeters
    config_t config;   // last accepted configuration of c and cdag
    weights_t w;       // weight values of the last accepted configuration
    u_frame_t u_frame; // frame of the last accepted configuraiton: just one time frame of a propagator

    qmc_config_data_t(atom_diag const &h_diag) : w{1., 1.} { u_frame = make_zero_propagator_frame(h_diag); }

    //qmc_config_data_t(params_t const &params, atom_diag const &h_diag, u_tau_t const &u_tau, block_gf_const_view<imtime> delta,
    //                  std::map<int, std::pair<int, int>> linindex);
    //int size() { return config.size(); }
  };

  struct hyb_adaptor_t {
    h_tau_t const &hyb_tau; // make a copy.
    std::map<int, std::pair<int, int>> const &linindex;

    hyb_adaptor_t(h_tau_t const &hyb_tau, std::map<int, std::pair<int, int>> const &linindex)
       : hyb_tau(std::move(hyb_tau)), linindex(std::move(linindex)) {}

    scalar_t operator()(double tau, int li, double tau_dag, int li_dag) const {
      auto [bl, in]         = linindex[li];
      auto [bl_dag, in_dag] = linindex[li_dag];

      EXPECTS(bl == bl_dag);
      double dtau = tau_dag - tau;
      if (dtau >= 0)
        return (hyb_tau[bl])(dtau)(in_dag, in);
      else
        return -(hyb_tau[bl])(hyb_tau[bl].domain().beta + dtau)(in_dag, in);
    }
  };

  struct qmc_params_t {

    // same for every inch step:
    hyb_adaptor_t hyb_adaptor;
    atom_diag const &h_diag; // Diagonalization of the atomic problem

    // updated at every inch step:
    u_tau_t const &u_tau; // THE propagator

    // different at every inch step:
    int inch_step;
    double tau_max;           // similar to beta, but configuration here does not always goes up to beta. 0 < tau_max <= beta
    double tau_split;         // in the inchworm, this should be the tau_max of the previous inching. 0 < tau_split <= tau_max
    bool use_bare_propagator; // true only for the first iteration of the inchworm calculation

    qmc_params_t(h_tau_t const &hyb_tau, std::map<int, std::pair<int, int>> const &linindex, atom_diag const &h_diag, u_tau_t const &u_tau,
                 int inch_step, double tau_max, double tau_split, bool use_bare_propagator)
       : hyb_adaptor(hyb_tau, linindex),
         h_diag(h_diag),
         u_tau(u_tau),
         inch_step(inch_step),
         tau_max(tau_max),
         tau_split(tau_split),
         use_bare_propagator(use_bare_propagator) {}
  };
} // namespace inchworm
