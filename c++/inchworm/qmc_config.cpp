#include "./qmc_config.hpp"

namespace inchworm {

  qmc_config_t::qmc_config_t(params_t const &params, atom_diag const &h_diag, block_gf_const_view<imtime> delta)
     : tau_seg(params.beta), h_diag(h_diag), delta(delta), last_accepted_U_frame(h_diag) {
    // ...
  }

} // namespace inchworm
