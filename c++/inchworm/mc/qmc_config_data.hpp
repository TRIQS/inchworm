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
    weights_t w;       // last accepted weight values
    u_frame_t u_frame; // last accepted frame: just one time frame of a propagator

    qmc_config_data_t(params_t const &params, atom_diag const &h_diag, u_tau_t const &u_tau, block_gf_const_view<imtime> delta,
                      std::map<int, std::pair<int, int>> linindex);
    //int size() { return config.size(); }
    //};

    //struct qmc_params_t {

    // same for every inch step:
    atom_diag const &h_diag; // Diagonalization of the atomic problem
                             //delta_adaptor_t delta_adaptor;

    block_gf<imtime> delta; // Hybridization function
    std::map<int, std::pair<int, int>>
       linindex; // structure that link the linear index from fops to the block/inner indices of the delta block function

    // updated at every inch step:
    u_tau_t const &u_tau; // THE propagator

    // changed at every inch step:
    int inch_step;
    double tau_max;           // similar to beta, but configuration here does not always goes up to beta. 0 < tau_max <= beta
    double tau_split;         // in the inchworm, this should be the tau_max of the previous inching. 0 < tau_split <= tau_max
    bool use_bare_propagator; // true only for the first iteration of the inchworm calculation
  };
} // namespace inchworm
