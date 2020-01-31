#include "./remove.hpp"
#include "../diagram/print.hpp"

namespace inchworm::moves {

  mc_weight_t remove::attempt() {

    //print_diag(data.get_time_diagram());
    int N        = data.size();
    bool success = false;
    while (not success) {
      int i     = rng(N);
      int i_dag = rng(N);
      success   = data.try_erase(i, i_dag);
      //std::printf("success? %d \n", success);
    }

    //std::printf("removing:");
    if (data.use_bare_propagator) {
      auto diagram = data.get_time_diagram();
      //print_configuration(diagram);
      new_w_hyb        = determinant(diagram);
      auto new_U_frame = propagator_product(data.h_diag, diagram, data.tau_max());
      new_w_loc        = new_U_frame.frobenius_norm();
    } else {
      std::vector<double> split_times = {data.tau_split()};
      auto diagram                    = data.get_time_diagram(split_times);
      new_w_hyb                       = inclusion_exclusion(diagram);
      auto new_U_frame                = propagator_product(data.U_tau, data.h_diag, diagram, data.tau_max());
      new_w_loc                       = new_U_frame.frobenius_norm();
    }

    auto w_hyb_ratio = new_w_hyb / data.last_accepted_w_hyb;
    auto w_loc_ratio = new_w_loc / data.last_accepted_w_loc;
    auto t_ratio     = std::pow(N / data.tau_max(), 2);

    //printf("helloaaaaaaaaaaiii02:  %f  \n", t_ratio * w_loc_ratio * w_hyb_ratio);
    if (false) {
      printf("new_w_hyb, new_w_loc, last_w_hyb, last_w_loc, t_ratio:  %f %f   %f %f   %f\n", new_w_hyb, new_w_loc, data.last_accepted_w_hyb,
             data.last_accepted_w_loc, t_ratio);
      fflush(stdout);
    }
    return t_ratio * w_loc_ratio * w_hyb_ratio;
  }

  mc_weight_t remove::accept() {
    std::printf("yes\n");
    print_configuration(data.get_time_diagram());
    data.last_accepted_w_hyb = new_w_hyb;
    data.last_accepted_w_loc = new_w_loc;
    data.update_accepted_lists();
    //data.last_accepted_U_frame = new_U_frame;

    return 1;
  }

  void remove::reject() {}

} // namespace inchworm::moves
