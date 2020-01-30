#include "./insert.hpp"

namespace inchworm::moves {

  mc_weight_t insert::attempt() {

    // Pick up the value of alpha and choose the operators
    int n_fops = qmc_config.h_diag.get_fops().size();
    int li1    = rng(n_fops);
    int li2    = rng(n_fops);

    // Choice of times for insertion. Find the time as double and them put them on the grid.
    auto tau1 = double(qmc_config.tau_seg.get_random_pt(rng));
    auto tau2 = double(qmc_config.tau_seg.get_random_pt(rng));

    int N        = qmc_config.size();
    bool success = false;
    while (not success) success = qmc_config.insert(tau1, li1, tau2, li2);

    if (qmc_config.use_bare_propagator) {
      new_w_hyb   = determinant(qmc_config.get_time_diagram());
      new_U_frame = propagator_product(qmc_config.h_diag, qmc_config.get_time_diagram(), qmc_config.tau_max());
    } else {
      std::vector<double> split_times = {qmc_config.tau_split()};
      new_w_hyb                       = inclusion_exclusion(qmc_config.get_time_diagram(split_times));
      new_U_frame                     = propagator_product(qmc_config.U_tau, qmc_config.h_diag, qmc_config.get_time_diagram(), qmc_config.tau_max());
    }

    auto w_hyb_ratio = new_w_hyb / qmc_config.last_accepted_w_hyb;

    // atomic weight
    new_w_loc        = new_U_frame.frobenius_norm();
    auto w_loc_ratio = new_w_loc / qmc_config.last_accepted_w_loc;

    // proposition probability
    auto t_ratio = std::pow(qmc_config.tau_max() / (N + 1), 2);
    return t_ratio * w_loc_ratio * w_hyb_ratio;
  }

  mc_weight_t insert::accept() {

    qmc_config.last_accepted_w_hyb   = new_w_hyb;
    qmc_config.last_accepted_w_loc   = new_w_loc;
    qmc_config.last_accepted_U_frame = new_U_frame;

    return 1;
  }

  void insert::reject() { qmc_config.erase_last(); }

} // namespace inchworm::moves
