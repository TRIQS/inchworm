#pragma once

#include "types.hpp"
#include "u_frame.hpp"

namespace inchworm {

  struct config_t {

    scalar_t imp_weight; // Impurity weight - Frobenius norm of the current (propagator or green function) frame
    scalar_t hyb_weight; // Hybridization weight - Value of the determinant in cthyb or its equivalent for the inchworm

    std::vector<fop_t> d_list, d_dag_list;                    // list of d/d_dag time ordered
    std::vector<std::vector<fop_t>> d_bl_list, d_dag_bl_list; // list of d/d_dag by block, insertion ordered

    int sign = 1; // The sign of the configuration

    // Create an empty configuration
    config_t(frame_t const &frame, gf_struct_t const &gf_struct)
       : imp_weight{frobenius_norm(frame)}, hyb_weight{1.0}, d_bl_list(gf_struct.size()), d_dag_bl_list(gf_struct.size()) {}

    long size() const { return d_list.size(); }
    long size(long bl) const { return d_bl_list[bl].size(); }

    bool try_insert(fop_t const &ddag, fop_t const &d);
    bool try_erase(long bl, long i_dag, long i);
    bool try_double_insert(fop_t const &d_dag1, fop_t const &d1, fop_t const &d_dag2, fop_t const &d2);
    bool try_double_erase(long bl1, long i1_dag, long i1, long bl2, long i2_dag, long i2);
  };

} // namespace inchworm
