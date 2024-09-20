#include "./impurity_product.hpp"
#include "./atom_diag.hpp"
#include "./u_frame.hpp"
#include "./interpolator.hpp"

namespace inchworm {

  bool has_zero_trace(atom_diag const &ad, diagram::time_diagram_t const &diagram) {
    if (diagram.size() == 0) return false;

    for (long initial_bl : range(ad.n_subspaces())) {
      long curr_bl = initial_bl;
      for (auto const &op : diagram.op_list) {
        curr_bl = (op.dag ? ad.cdag_connection(op.linear_index, curr_bl) : ad.c_connection(op.linear_index, curr_bl));
        if (curr_bl == -1) break;
      }
      // Any non-void block gives us a finite trace contribution
      if (curr_bl != -1) return false;
    }
    return true;
  }

  u_partial_t impurity_product(atom_diag const &ad, diagram::time_diagram_t const &diagram, double tau_max, double tau_min, double energy_shift,
                               interpolator_t<scalar_t> const *const u_interpolator_p) {
    EXPECTS(tau_max > tau_min);

    // Filter out all operators in the time-window [tau_min, tau_max]
    std::vector<fop_t> op_lst;
    op_lst.reserve(diagram.op_list.size());
    for (auto const &op : diagram.op_list) {
      if ((op.tau >= tau_min) and (op.tau <= tau_max)) op_lst.push_back(op);
    }

    // Helper to calculate u_tau for both inchworm and cthyb case
    auto u_tau = [&](int bl, double tau) -> matrix_t {
      if (u_interpolator_p) { // inchworm case, interpolate
        return (*u_interpolator_p)(bl, tau);
      } else { // cthyb case
        auto bl_size = ad.get_subspace_dim(bl);
        auto res     = matrix_t::zeros({bl_size, bl_size});
        for (auto j : range(bl_size)) res(j, j) = std::exp(-tau * (ad.get_eigenvalue(bl, j) + energy_shift));
        return res;
      }
    };

    // Calculate full operator product
    u_partial_t u_partial(ad.n_subspaces());
    if (op_lst.empty()) {
      for (long initial_bl : range(ad.n_subspaces())) u_partial[initial_bl] = {initial_bl, u_tau(initial_bl, tau_max - tau_min)};
      return u_partial;
    }

    for (long initial_bl : range(ad.n_subspaces())) {

      // Short-circuit if the product contains a void block, i.e. -1
      long curr_bl = initial_bl;
      for (auto const &op : op_lst) {
        curr_bl = (op.dag ? ad.cdag_connection(op.linear_index, curr_bl) : ad.c_connection(op.linear_index, curr_bl));
        if (curr_bl == -1) break;
      }
      if (curr_bl == -1) {
        u_partial[initial_bl] = {-1, {}};
        continue; // Skip matrix product
      }

      // Perform Matrix multiplication
      // u(tau_max - tau_{n-1}) op_{n-1} u(tau_{n-1}-tau_{n-2}) op_{n-2} ... op_1 u(tau_1 - tau_0) op_0 u(tau_0 - tau_min)
      auto mat_bl = u_tau(initial_bl, op_lst.front().tau - tau_min);
      curr_bl     = initial_bl;
      for (auto [i, op] : enumerate(op_lst)) {
        mat_bl  = (op.dag ? ad.cdag_matrix(op.linear_index, curr_bl) : ad.c_matrix(op.linear_index, curr_bl)) * mat_bl;
        curr_bl = (op.dag ? ad.cdag_connection(op.linear_index, curr_bl) : ad.c_connection(op.linear_index, curr_bl));

        if (i == op_lst.size() - 1)
          mat_bl = u_tau(curr_bl, tau_max - op.tau) * mat_bl;
        else
          mat_bl = u_tau(curr_bl, op_lst[i + 1].tau - op.tau) * mat_bl;
      }
      u_partial[initial_bl] = {curr_bl, mat_bl};
    }
    return u_partial;
  }

} // namespace inchworm
