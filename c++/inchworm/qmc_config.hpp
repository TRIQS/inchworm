#pragma once
#include "./types.hpp"
#include "./params.hpp"
#include "./configuration.hpp"
//#include "./diagram/diagram.hpp"
//#include "./impurity_product.hpp"
#include <triqs/utility/time_pt.hpp>
#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/utility/serialization.hpp>
#include <triqs/det_manip.hpp>

namespace inchworm {
  using triqs::utility::time_segment;

  /// The Monte-Carlo Configuration Class
  struct qmc_config_t {

    configuration config;                        // Configuration
    time_segment tau_seg;                        // discretized time segment
    std::map<std::pair<int, int>, int> linindex; // Linear index constructed from block and inner indices
    atom_diag const &h_diag;                     // Diagonalization of the atomic problem

    std::vector<int> n_inner;
    block_gf<imtime> delta; // Hybridization function

    /// This callable object adapts the Delta function for the call of the det.
    struct delta_block_adaptor {
      gf<imtime, delta_target_t> delta_block; // make a copy. Needed in the real case anyway.

      delta_block_adaptor(gf<imtime, delta_target_t> delta_block) : delta_block(std::move(delta_block)) {}
      delta_block_adaptor(delta_block_adaptor const &) = default;
      delta_block_adaptor(delta_block_adaptor &&)      = default;
      delta_block_adaptor &operator=(delta_block_adaptor const &) = delete;
      delta_block_adaptor &operator=(delta_block_adaptor &&) = default;

      det_scalar_t operator()(std::pair<time_pt, int> const &x, std::pair<time_pt, int> const &y) const {
        det_scalar_t res = delta_block[closest_mesh_pt(double(x.first - y.first))](x.second, y.second);
        return (x.first >= y.first ? res : -res); // x,y first are time_pt, wrapping is automatic in the - operation, but need to
                                                  // compute the sign
      }

      friend void swap(delta_block_adaptor &dba1, delta_block_adaptor &dba2) noexcept { swap(dba1.delta_block, dba2.delta_block); }
    };

    std::vector<triqs::det_manip::det_manip<delta_block_adaptor>> dets; // The determinants
    int current_sign, old_sign;                                         // Permutation prefactor
    h_scalar_t atomic_weight;                                           // The current value of the trace or norm
    h_scalar_t atomic_reweighting;                                      // The current value of the reweighting

    /*
    time_diagram_t diagram;                           // Configuration
    //mutable impurity_product imp_prod;                // Calculator of THE product
    //block_gf<imtime, delta_target_t> delta; // Hybridization function

    //Constructor
*/
    //qmc_config_t(params_t const &params);
    qmc_config_t(params_t const &params, atom_diag const &h_diag, std::map<std::pair<int, int>, int> linindex, block_gf_const_view<imtime> delta,
                 std::vector<int> n_inner);
  };
} // namespace inchworm
