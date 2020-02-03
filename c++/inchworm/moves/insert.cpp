#include "./insert.hpp"
#include "../diagram/print.hpp"

namespace inchworm::moves {

  mc_weight_t insert::attempt() {

    //print_diag(data.get_time_diagram());
    int n_fops = (data.h_diag.get_fops()).size();
    int li1    = rng(n_fops);
    int li2    = rng(n_fops);

    //bool success = false;
    //while (not success) {
    // Choice of times for insertion. Try until success.
    double tau1 = rng(data.tau_max());
    double tau2 = rng(data.tau_max());
    if (not data.try_insert(tau1, li1, tau2, li2)) return 0; //data is not modified in this case

    //std::printf("success? %d \n", success);
    //}
    int N = data.size();

    //std::printf("\ninserting:");
    if (data.use_bare_propagator) {
      auto diagram = diagram::time_diagram_t{config.c_list, config.cdag_list, {}}; // make a free function (not member of data)
      //print_configuration(diagram);
      new_w_hyb        = determinant(diagram);
      auto new_U_frame = propagator_product(data.h_diag, diagram, data.tau_max());
    } else {
      std::vector<double> split_times = {data.tau_split()}; //put in constructor of move
      //auto diagram                    = data.get_time_diagram(split_times); // make a free function (not member of data)

      auto diagram     = diagram::time_diagram_t{config.c_list, config.cdag_list, split_times}; // make a free function (not member of data)
      new_w_hyb        = diagram::inclusion_exclusion(diagram);
      auto new_U_frame = propagator_product(data.U_tau, data.h_diag, diagram, data.tau_max());
    }
    new_w_loc = new_U_frame.frobenius_norm();

    auto w_hyb_ratio = new_w_hyb / data.last_accepted_w_hyb;
    auto w_loc_ratio = new_w_loc / data.last_accepted_w_loc;
    auto t_ratio     = std::pow(data.tau_max() / (N + 1), 2);

    //printf("helloaaaaaaaaaaiii02:  %f  \n", t_ratio * w_loc_ratio * w_hyb_ratio);
    if (false) {
      printf("new_w_hyb, new_w_loc, last_w_hyb, last_w_loc, t_ratio:  %f %f   %f %f   %f\n", new_w_hyb, new_w_loc, data.last_accepted_w_hyb,
             data.last_accepted_w_loc, t_ratio);
      fflush(stdout);
    }
    return t_ratio * w_loc_ratio * w_hyb_ratio;
  }

  //
  mc_weight_t insert::accept() {
    //std::printf("yes\n");
    print_configuration(data.get_time_diagram());
    data.last_accepted_w_hyb = new_w_hyb;
    data.last_accepted_w_loc = new_w_loc;
    data.update_accepted_lists();

    return 1;
  }

} // namespace inchworm::moves
