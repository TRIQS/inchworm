# Generated automatically using the command :
# c++2py ../../c++/inchworm/atom_diag.hpp -p --members_read_only -N inchworm -a inchworm -m ad_tools -o ad_tools --moduledoc="Inchworm atom_diag tools" -C triqs -C nda --cxxflags="-std=c++20" --target_file_only
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "ad_tools", doc = r"Inchworm atom_diag tools", app_name = "inchworm")

# Imports
module.add_imports(*['triqs.atom_diag', 'triqs.gf', 'triqs.operators'])

# Add here all includes
module.add_include("inchworm/atom_diag.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/pair.hpp>
#include <cpp2py/converters/std_array.hpp>
#include <cpp2py/converters/string.hpp>
#include <cpp2py/converters/vector.hpp>
#include <nda_py/cpp2py_converters.hpp>
#include <triqs/cpp2py_converters/gf.hpp>
#include <triqs/cpp2py_converters/operators_real_complex.hpp>
#include <triqs/cpp2py_converters/real_or_complex.hpp>

using namespace inchworm;
""")


module.add_function ("triqs::operators::many_body_operator inchworm::create_effective_hyb (triqs::hilbert_space::gf_struct_t gf_struct)", doc = r"""""")

module.add_function ("inchworm::u_partial_t inchworm::get_op_block_matrix (inchworm::atom_diag ad, std::string bl_name, int idx, bool op_dag)", doc = r"""""")

module.add_function ("inchworm::frame_t inchworm::make_g_frame_from_l_and_r (inchworm::atom_diag ad_imp, triqs::hilbert_space::gf_struct_t gf_struct, inchworm::u_partial_t l, inchworm::u_partial_t r)", doc = r"""""")

module.add_function ("inchworm::frame_t inchworm::make_bare_u_frame (inchworm::atom_diag ad, double tau, bool set_gs_to_0 = false)", doc = r"""""")

module.add_function ("inchworm::frame_t inchworm::make_bare_g_frame (inchworm::atom_diag ad_imp, inchworm::u_tau_t u_tau, triqs::hilbert_space::gf_struct_t gf_struct, double tau_split, double beta)", doc = r"""""")

module.add_function ("inchworm::frame_t inchworm::partial_trace_bath (inchworm::atom_diag ad_tot, inchworm::atom_diag ad_target, inchworm::atom_diag ad_bath, double beta, double tau)", doc = r"""Calculate U(tau) = Tr_bath [ exp[-(beta - tau) * H_bath] * exp(-tau * H_tot) ] / Z_bath

 Within atom_diag, perform partial trace over indices above nfops_imp, i.e. the bath degrees of freedom.
 Only the nfops_imp first degrees of freedom will be preserved.

Parameters
----------
ad_tot
     atom_diag of the full Hamiltonian H = H_loc + H_bath + H_hyb.

ad_imp
     atom_diag of the local Hamiltonian H_loc.

ad_bath
     atom_diag of the bath Hamiltonian H_bath (in the full basis).

beta
     inverse temperature.

tau
     0 < tau < beta.

Returns
-------
out
     The partial sum matrix of a function H. The result is a block diagonal matrix, with blocks and indices in the same order as in ad_imp.""")

module.add_function ("inchworm::u_tau_t inchworm::make_ED_propagator (inchworm::atom_diag ad_tot, inchworm::atom_diag ad_imp, inchworm::atom_diag ad_bath, double beta, int n_tau)", doc = r"""""")



module.generate_code()