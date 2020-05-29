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

  u_tau_t make_propagator(atom_diag const &h_diag, double beta, int n_tau);
  u_tau_t make_ED_propagator(atom_diag const &ad_tot, atom_diag const &ad_atom, atom_diag const &ad_bath, double beta, int n_tau);

  //
  u_frame_t make_zeroth_order(atom_diag const &ad, double tau, double tau_split = 0.0, u_tau_t const *const u_tau_p = nullptr);

  /// function that calculate the product: u_frame = u(tau_0) op u(tau_1-tau_0) op u(tau_2-tau_1) op u(tau_3-tau_2) ... op u(tau-tau_n)
  /// where op is either c_dag or c operator, depending on the configuration
  /** 
   * @param ad atom_diag of the system considered here.
   * @param diagram configuration of the n operators (op) of the present monte carlo step.
   * @param tau_max time of the u_frame_t calculated here.
   * @param tau_split time interval of the previously calculated propagator (u_tau).
   * @param u_tau_p pointer to the full propagator (u_tau) calculated up until this point (0 < tau < tau_split). 
   * @return u_frame_t, at time tau_max, resulting from this product.
   */
  u_frame_t propagator_product(atom_diag const &ad, time_diagram_t const &diagram, double tau_max, double tau_split = 0.0,
                               u_tau_t const *const u_tau_p = nullptr);

  inline std::ostream &operator<<(std::ostream &out, u_frame_t const &u_frame) {
    out << "propagator_frame (size: " << u_frame.size() << ")\n";
    for (int bl = 0; bl < u_frame.size(); bl++) { out << u_frame[bl] << "\n"; }
    return out;
  }

  constexpr int MAX_ORDER = 7; // just for printing
  struct single_step_results_t {
    double average_k = 0.0;
    u_frame_t u_frame;
    u_frame_t u_frame_0th_order;
    std::vector<double> u_expansion_order;
    std::vector<int> samples_expansion_order;
    single_step_results_t(atom_diag const &h_diag) : u_expansion_order(MAX_ORDER, 0), samples_expansion_order(MAX_ORDER, 0) {
      u_frame           = make_zero_propagator_frame(h_diag);
      u_frame_0th_order = u_frame;
    };

    void normalize(double normalization_cte) {
      for (auto &Bl : u_frame) Bl /= normalization_cte;
      for (auto &Bl : u_frame_0th_order) Bl /= normalization_cte;
      for (auto &o : u_expansion_order) o /= normalization_cte;
    };

    void print();
   
    
  };

} // namespace inchworm
