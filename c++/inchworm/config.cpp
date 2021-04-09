#include "./config.hpp"

namespace inchworm {

  inline void insert_sorted(auto &vec, auto const &el) { vec.insert(std::upper_bound(vec.begin(), vec.end(), el), el); }

  bool config_t::try_insert(fop_t const &d_dag, fop_t const &d) {
    EXPECTS(d_dag.bl == d.bl);

    for (int i = 0; i < size() - 1; i++)
      if ((d_list[i].tau == d.tau) or (d_dag_list[i].tau == d_dag.tau)) return false;

    insert_sorted(d_list, d);
    insert_sorted(d_dag_list, d_dag);

    d_bl_list[d.bl].push_back(d);
    d_dag_bl_list[d_dag.bl].push_back(d_dag);

    return true;
  }

  bool config_t::try_erase(long bl, long i_dag, long i) {
    EXPECTS(i < size(bl) && i_dag < size(bl));

    std::erase(d_list, d_bl_list[bl][i]);
    std::erase(d_dag_list, d_dag_bl_list[bl][i_dag]);

    d_bl_list[bl].erase(d_bl_list[bl].begin() + i);
    d_dag_bl_list[bl].erase(d_dag_bl_list[bl].begin() + i_dag);

    return true;
  }

  bool config_t::try_double_insert(fop_t const &d_dag1, fop_t const &d1, fop_t const &d_dag2, fop_t const &d2) {
    EXPECTS(d_dag1.bl == d1.bl && d_dag2.bl == d2.bl);

    // Make sure that no two times are equal
    // FIXME Protect equal time operator insertions
    for (int i = 0; i < size() - 1; i++)
      if ((d_list[i].tau == d1.tau) or (d_dag_list[i].tau == d_dag1.tau) or (d_list[i].tau == d2.tau) or (d_dag_list[i].tau == d_dag2.tau)
          or (d1.tau == d2.tau) or (d_dag1.tau == d_dag2.tau))
        return false;

    insert_sorted(d_list, d1);
    insert_sorted(d_list, d2);
    insert_sorted(d_dag_list, d_dag1);
    insert_sorted(d_dag_list, d_dag2);

    d_bl_list[d1.bl].push_back(d1);
    d_bl_list[d2.bl].push_back(d2);
    d_dag_bl_list[d_dag1.bl].push_back(d_dag1);
    d_dag_bl_list[d_dag2.bl].push_back(d_dag2);

    return true;
  }

  bool config_t::try_double_erase(long bl1, long i1_dag, long i1, long bl2, long i2_dag, long i2) {
    EXPECTS(i1 < size(bl1) and i1_dag < size(bl1));
    EXPECTS(i2 < size(bl2) and i2_dag < size(bl2));

    if (bl1 == bl2 and (i1 == i2 or i1_dag == i2_dag)) return false;

    std::erase(d_list, d_bl_list[bl1][i1]);
    std::erase(d_list, d_bl_list[bl2][i2]);
    std::erase(d_dag_list, d_dag_bl_list[bl1][i1_dag]);
    std::erase(d_dag_list, d_dag_bl_list[bl2][i2_dag]);

    if (i1 < i2 && bl1 == bl2) std::swap(i1, i2);
    d_bl_list[bl1].erase(d_bl_list[bl1].begin() + i1);
    d_bl_list[bl2].erase(d_bl_list[bl2].begin() + i2);

    if (i1_dag < i2_dag && bl1 == bl2) std::swap(i1_dag, i2_dag);
    d_dag_bl_list[bl1].erase(d_dag_bl_list[bl1].begin() + i1_dag);
    d_dag_bl_list[bl2].erase(d_dag_bl_list[bl2].begin() + i2_dag);

    return true;
  }

} // namespace inchworm
