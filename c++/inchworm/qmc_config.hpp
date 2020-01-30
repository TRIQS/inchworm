#pragma once
#include "./types.hpp"
#include "./params.hpp"
#include "./util.hpp"
#include "./diagram/diagram.hpp"
#include "./impurity_product.hpp"
#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/atom_diag/functions.hpp>
#include <triqs/utility/serialization.hpp>
#include <triqs/det_manip.hpp>

namespace inchworm {
  using triqs::utility::time_segment;

  /// The Monte-Carlo Configuration Class
  struct qmc_config_t {

    std::vector<time_and_index_t> c_list, cdag_list; // list of c/cdag not time ordered, but different
    double tau_max_;                                 // similar to beta, but configuration here does not always goes up to beta. 0 < tau_max_ <= beta
    double tau_split_;                               // should be the tau_max_ of the previous inching.  0 < tau_split_ <= tau_max_

    //configuration config;    // Configuration
    atom_diag const &h_diag; // Diagonalization of the atomic problem
    block_gf<imtime> delta;  // Hybridization function

    int current_sign, old_sign;     // Permutation prefactor
    h_scalar_t last_accepted_w_loc; // atomic weight (Frobenius norm of the current propagator frame)
    h_scalar_t last_accepted_w_hyb; // value of the determinant in cthyb or its equivalent for the inchworm
    propagator_frame last_accepted_U_frame;
    u_tau_t U_tau; // THE propagator

    bool use_bare_propagator; //true only for the first iteration of the inchworm calculation

    //mutable impurity_product imp_prod;                // Calculator of THE product

    qmc_config_t(params_t const &params, atom_diag const &h_diag, block_gf_const_view<imtime> delta);

    double tau_max();
    double tau_split();
    int size();

    bool insert(double tau, int linear_index, double tau_dag, int linear_index_dag);
    bool erase(int i, int i_dag);
    bool erase_last();
    void clear();

    time_diagram_t get_time_diagram(std::vector<double> const &split_times);

    time_diagram_t get_time_diagram();
  };
} // namespace inchworm
