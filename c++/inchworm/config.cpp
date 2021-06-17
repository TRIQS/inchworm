#include "./config.hpp"

namespace inchworm {

  inline void insert_sorted(auto &vec, auto const &el) { vec.insert(std::upper_bound(vec.begin(), vec.end(), el), el); }

  bool config_t::try_insert(fop_t const &d_dag, fop_t const &d) {
    EXPECTS(d_dag.bl == d.bl);

    if (std::count(split_times.cbegin(), split_times.cend(), d.tau) > 0) return false;
    if (std::count(split_times.cbegin(), split_times.cend(), d_dag.tau) > 0) return false;

    split_times.push_back(d.tau);
    split_times.push_back(d_dag.tau);

    insert_sorted(d_list, d);
    insert_sorted(d_dag_list, d_dag);

    d_bl_list[d.bl].push_back(d);
    d_dag_bl_list[d_dag.bl].push_back(d_dag);

    return true;
  }

  bool config_t::try_erase(long bl, long i_dag, long i) {
    EXPECTS(i < size(bl) && i_dag < size(bl));

    std::erase(split_times, d_bl_list[bl][i].tau);
    std::erase(split_times, d_dag_bl_list[bl][i_dag].tau);

    std::erase(d_list, d_bl_list[bl][i]);
    std::erase(d_dag_list, d_dag_bl_list[bl][i_dag]);

    d_bl_list[bl].erase(d_bl_list[bl].begin() + i);
    d_dag_bl_list[bl].erase(d_dag_bl_list[bl].begin() + i_dag);

    return true;
  }

} // namespace inchworm
