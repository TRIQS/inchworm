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
    if ((size() <= i) or (size() <= i_dag)) {
      std::printf("heille.\n");
      return false;
    }
    d_list.erase(d_list.begin() + i);
    d_dag_list.erase(d_dag_list.begin() + i_dag);
    return true;
  }

  bool config_t::try_double_insert(double tau1, int linear_index1, double tau1_dag, int linear_index1_dag, double tau2, int linear_index2,
                                   double tau2_dag, int linear_index2_dag) {
    for (int i = 0; i < size() - 1; i++)
      if ((d_list[i].tau == tau1) or (d_dag_list[i].tau == tau1_dag) or (d_list[i].tau == tau2) or (d_dag_list[i].tau == tau2_dag) or (tau1 == tau2)
          or (tau1_dag == tau2_dag))
        return false;
    d_list.push_back({tau1, linear_index1});
    d_list.push_back({tau2, linear_index2});
    d_dag_list.push_back({tau1_dag, linear_index1_dag});
    d_dag_list.push_back({tau2_dag, linear_index2_dag});
    return true;
  }

  bool config_t::try_double_erase(int i, int i_dag, int j, int j_dag) {
    //std::printf("try_erase? %d   %d %d \n", size(), i, i_dag);
    if (i == j or i_dag == j_dag) return false;
    if ((size() <= i) or (size() <= i_dag) or (size() <= j_dag) or (size() <= j_dag)) {
      std::printf("heille.\n");
      return false;
    } // FIXME: refactoring? (expects)
    if (i < j) std::swap(i, j);
    d_list.erase(d_list.begin() + i);
    d_list.erase(d_list.begin() + j);

    if (i_dag < j_dag) std::swap(i_dag, j_dag);
    d_dag_list.erase(d_dag_list.begin() + i_dag);
    d_dag_list.erase(d_dag_list.begin() + j_dag);
    return true;
  }

} // namespace inchworm
