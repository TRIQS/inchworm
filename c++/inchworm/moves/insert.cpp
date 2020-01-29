#include "./update.hpp"

namespace inchworm::moves {

  mc_weight_t insert::attempt() {

    // Pick up the value of alpha and choose the operators
    //auto rs1 = rng(block_size), rs2 = rng(block_size);

    int bl1 = rng(qmc_config.n_inner.size());
    int bl2 = rng(qmc_config.n_inner.size()); //at some point, should be equal to bl1.

    int rs1 = rng(qmc_config.n_inner[bl1]);
    int rs2 = rng(qmc_config.n_inner[bl2]);

    // Choice of times for insertion. Find the time as double and them put them on the grid.
    auto tau1 = qmc_config.tau_seg.get_random_pt(rng);
    auto tau2 = qmc_config.tau_seg.get_random_pt(rng);
    // FIXME need to check if any tau is already included in the config.

    // record the length of the proposed insertion
    auto dtau = double(tau2 - tau1);
    //if (histo_proposed) *histo_proposed << dtau;

    // Computation of det ratio
    //auto det_ratio = det.try_insert(num_c_dag, num_c, {tau1, op1.inner_index}, {tau2, op2.inner_index});

    // proposition probability
    //mc_weight_t t_ratio = std::pow(block_size * config.beta() / double(det.size() + 1), 2);

    return 0;
  }

  mc_weight_t insert::accept() {
    /*
    // insert in the tree
    data.imp_trace.confirm_insert();

    // insert in the configuration
    config.insert(tau1, op1);
    config.insert(tau2, op2);
    config.finalize();

    // insert in the determinant
    data.dets[block_index].complete_operation();
    data.update_sign();
    data.atomic_weight      = new_atomic_weight;
    data.atomic_reweighting = new_atomic_reweighting;
    if (histo_accepted) *histo_accepted << dtau;

#ifdef EXT_DEBUG
    std::cerr << "* Move move_insert_c_cdag accepted" << std::endl;
    std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    check_det_sequence(data.dets[block_index], config.get_id());
#endif

    return data.current_sign / data.old_sign;
*/
    return {};
  }

  void insert::reject() {
    /*
    config.finalize();
    data.imp_trace.cancel_insert();
    data.dets[block_index].reject_last_try();

#ifdef EXT_DEBUG
    std::cerr << "* Move move_insert_c_cdag rejected" << std::endl;
    std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    check_det_sequence(data.dets[block_index], config.get_id());
#endif
*/
  }

} // namespace inchworm::moves
