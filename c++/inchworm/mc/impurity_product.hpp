#pragma once
#include "../u_frame.hpp"
#include "../util.hpp"
#include "../types.hpp"
#include "../params.hpp"
#include "../diagram/diagram.hpp"

//#include <inchworm/solver_core.hpp>
#include <triqs/gfs.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/utility/macros.hpp>
#include <triqs/arrays.hpp>

namespace inchworm {
  using time_diagram_t = diagram::time_diagram_t;

  // FIXME: these three functions should be put in a seperated file

  // Make an empty propagator (green function) with the same structure as the one in atom_diag:
  u_tau_t make_propagator(atom_diag const &ad_imp, double beta, int n_tau);

  // Make an exact diagonalization propagator U = Trace_B [exp(-H_bath *(beta-tau)) exp(-H_tot*tau)  ]  /  Trace_B [ exp(-H_bath*beta) ]
  u_tau_t make_ED_propagator(atom_diag const &ad_tot, atom_diag const &ad_imp, atom_diag const &ad_bath, double beta, int n_tau);

  // zeroth order is :
  // U_0(tau) = exp(-H_imp*tau) for cthyb
  // U_0(tau) = U(tau-tau_split) U(tau_split) for inchworm
  //
  frame_t make_zeroth_order_frame(atom_diag const &ad, double tau, double tau_split = 0.0, u_tau_t const *const u_tau_p = nullptr);

  /// function that calculate the product: u_frame = u(tau_0) op u(tau_1-tau_0) op u(tau_2-tau_1) op u(tau_3-tau_2) ... op u(tau-tau_n)
  /// where op is either c_dag or c operator, depending on the configuration
  /** 
   * @param ad atom_diag of the system considered here.
   * @param diagram configuration of the n operators (op) of the present monte carlo step.
   * @param tau_max time of the frame_t calculated here.
   * @param tau_split time interval of the previously calculated propagator (u_tau).
   * @param u_tau_p pointer to the full propagator (u_tau) calculated up until this point (0 < tau < tau_split). 
   * @return frame_t, at time tau_max, resulting from this product.
   */
  u_partial_t impurity_product(atom_diag const &ad, time_diagram_t const &diagram, double tau0, double tau1, u_tau_t const *const u_tau_p = nullptr);

  inline std::vector<matrix_t> make_frame(std::vector<long> shape_of_frame) {
    std::vector<matrix_t> res;
    for (auto n : shape_of_frame) { res.push_back(matrix_t(n, n)); }
    return res;
  }

  constexpr int MAX_ORDER = 7; // just for printing

  // structure to gather result of one Monte Carlo run:
  struct single_step_results_t {
    double average_k = 0.0;
    frame_t frame;
    frame_t frame_0th_order;
    std::vector<double> expansion_order;
    std::vector<int> samples_expansion_order;
    single_step_results_t(std::vector<long> shape_of_frame) : expansion_order(MAX_ORDER, 0), samples_expansion_order(MAX_ORDER, 0) {
      frame           = make_frame(shape_of_frame);
      frame_0th_order = frame;
    };

    void normalize(double normalization_cte) {
      for (auto &Bl : frame) Bl /= normalization_cte;
      for (auto &Bl : frame_0th_order) Bl /= normalization_cte;
      for (auto &o : expansion_order) o /= normalization_cte;
    };

    void print();
  };

} // namespace inchworm
