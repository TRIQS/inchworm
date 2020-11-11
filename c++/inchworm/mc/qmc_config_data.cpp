#include "./qmc_config_data.hpp"
#include "impurity_product.hpp"

namespace inchworm {

  bool config_t::try_insert(fop_t const &d_dag, fop_t const &d) {
    for (int i = 0; i < size() - 1; i++)
      if ((d_list[i].tau == d.tau) or (d_dag_list[i].tau == d_dag.tau)) return false;
    d_list.push_back(d);
    d_dag_list.push_back(d_dag);
    return true;
  }

  bool config_t::try_erase(int i, int i_dag) {
    if ((size() <= i) or (size() <= i_dag)) {
      std::printf("heille.\n");
      return false;
    }
    d_list.erase(d_list.begin() + i);
    d_dag_list.erase(d_dag_list.begin() + i_dag);
    return true;
  }

  bool config_t::try_double_insert(fop_t const &d_dag1, fop_t const &d1, fop_t const &d_dag2, fop_t const &d2) {
    for (int i = 0; i < size() - 1; i++)
      if ((d_list[i].tau == d1.tau) or (d_dag_list[i].tau == d_dag1.tau) or (d_list[i].tau == d2.tau) or (d_dag_list[i].tau == d_dag2.tau)
          or (d1.tau == d2.tau) or (d_dag1.tau == d_dag2.tau))
        return false;
    d_list.push_back(d1);
    d_list.push_back(d2);
    d_dag_list.push_back(d_dag1);
    d_dag_list.push_back(d_dag2);
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
