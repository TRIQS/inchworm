#include "./qmc_config.hpp"

namespace inchworm {

  qmc_config_t::qmc_config_t(params_t const &params, atom_diag const &h_diag, std::map<std::pair<int, int>, int> linindex,
                             block_gf_const_view<imtime> delta, std::vector<int> n_inner)
     : config(params.beta), tau_seg(params.beta), h_diag(h_diag), n_inner(n_inner) {
    // ...
  }

} // namespace inchworm
