#include "./qmc_config_data.hpp"
#include "impurity_product.hpp"

namespace inchworm {

  bool config_t::try_insert(double tau, int linear_index, double tau_dag, int linear_index_dag) {
    //std::printf("try_erase? %d   %2.4f %d  %2.4f %d \n", size(),  tau, linear_index, tau_dag, linear_index_dag);
    for (int i = 0; i < size() - 1; i++)
      if ((d_list[i].tau == tau) or (d_dag_list[i].tau == tau_dag)) return false;
    d_list.push_back({tau, linear_index}); // FIXME upper_bound ordering
    d_dag_list.push_back({tau_dag, linear_index_dag});
    return true;
  }

  bool config_t::try_erase(int i, int i_dag) {
    //std::printf("try_erase? %d   %d %d \n", size(), i, i_dag);
    if ((size() <= i) or (size() <= i_dag)){ std::printf("heille.\n"); return false;}
    d_list.erase(d_list.begin() + i);
    d_dag_list.erase(d_dag_list.begin() + i_dag);
    return true;
  }

} // namespace inchworm
