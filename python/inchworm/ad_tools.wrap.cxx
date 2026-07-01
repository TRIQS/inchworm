
// C.f. https://numpy.org/doc/1.21/reference/c-api/array.html#importing-the-api
#define PY_ARRAY_UNIQUE_SYMBOL _cpp2py_ARRAY_API
#ifndef CLAIR_C2PY_WRAP_GEN
#ifdef __clang__
// #pragma clang diagnostic ignored "-W#warnings"
#endif
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wcast-function-type"
#pragma GCC diagnostic ignored "-Wcpp"
#endif

#define C2PY_VERSION_MAJOR 0
#define C2PY_VERSION_MINOR 1

#include <c2py/c2py.hpp>

using c2py::operator""_a;

// ==================== enums =====================

// ==================== module classes =====================

// ==================== module functions ====================

// create_effective_hyb
static auto const _c2py_fun_0 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](const triqs::gfs::gf_struct_t &gf_struct) { return inchworm::create_effective_hyb(gf_struct); }, "gf_struct")};

// get_op_block_matrix
static auto const _c2py_fun_1 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](const inchworm::atom_diag &ad, const std::string &bl_name, long idx,
                                         bool op_dag) { return inchworm::get_op_block_matrix(ad, bl_name, idx, op_dag); },
                                      "ad", "bl_name", "idx", "op_dag")};

// make_ED_propagator
static auto const _c2py_fun_2 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](const inchworm::atom_diag &ad_tot, const inchworm::atom_diag &ad_imp, const inchworm::atom_diag &ad_bath,
                                         double beta, long n_tau) { return inchworm::make_ED_propagator(ad_tot, ad_imp, ad_bath, beta, n_tau); },
                                      "ad_tot", "ad_imp", "ad_bath", "beta", "n_tau")};

// make_bare_g_frame
static auto const _c2py_fun_3 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](const inchworm::atom_diag &ad_imp, const inchworm::u_tau_t &u_tau, const triqs::gfs::gf_struct_t &gf_struct, double tau_split,
                 double beta) { return inchworm::make_bare_g_frame(ad_imp, u_tau, gf_struct, tau_split, beta); },
              "ad_imp", "u_tau", "gf_struct", "tau_split", "beta")};

// make_bare_u_frame
static auto const _c2py_fun_4 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](const inchworm::atom_diag &ad, double tau) { return inchworm::make_bare_u_frame(ad, tau); }, "ad", "tau")};

// make_g_frame_from_l_and_r
static auto const _c2py_fun_5 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](const inchworm::atom_diag &ad_imp, const triqs::gfs::gf_struct_t &gf_struct, const inchworm::u_partial_t &l,
                                         const inchworm::u_partial_t &r) { return inchworm::make_g_frame_from_l_and_r(ad_imp, gf_struct, l, r); },
                                      "ad_imp", "gf_struct", "l", "r")};

// partial_trace_bath
static auto const _c2py_fun_6 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](const inchworm::atom_diag &ad_tot, const inchworm::atom_diag &ad_target, const inchworm::atom_diag &ad_bath,
                                         double beta, double tau) { return inchworm::partial_trace_bath(ad_tot, ad_target, ad_bath, beta, tau); },
                                      "ad_tot", "ad_target", "ad_bath", "beta", "tau")};

static const auto _c2py_doc_0 = _c2py_fun_0.doc(R"DOC()DOC");
static const auto _c2py_doc_1 = _c2py_fun_1.doc(R"DOC()DOC");
static const auto _c2py_doc_2 = _c2py_fun_2.doc(R"DOC()DOC");
static const auto _c2py_doc_3 = _c2py_fun_3.doc(R"DOC()DOC");
static const auto _c2py_doc_4 = _c2py_fun_4.doc(R"DOC()DOC");
static const auto _c2py_doc_5 = _c2py_fun_5.doc(R"DOC()DOC");
static const auto _c2py_doc_6 = _c2py_fun_6.doc(R"DOC(
Calculate U(tau) = Tr_bath [ exp[-(beta - tau) * H_bath] * exp(-tau * H_tot) ] / Z_bath

Within atom_diag, perform partial trace over indices above nfops_imp, i.e. the bath degrees of freedom.
Only the nfops_imp first degrees of freedom will be preserved.

Parameters
----------
ad_tot : {par_0}
   atom_diag of the full Hamiltonian H = H_loc + H_bath + H_hyb.
ad_imp : {par_1}
   atom_diag of the local Hamiltonian H_loc.
ad_bath : {par_2}
   atom_diag of the bath Hamiltonian H_bath (in the full basis).
beta : {par_3}
   inverse temperature.
tau : {par_4}
   0 < tau < beta.

Returns
-------
{ret_0}
   The partial sum matrix of a function H. The result is a block diagonal matrix, with blocks and indices in the same order as in ad_imp.
)DOC",
                                                {{c2py::python_typename<const inchworm::atom_diag &>()},
                                                 {},
                                                 {c2py::python_typename<const inchworm::atom_diag &>()},
                                                 {c2py::python_typename<double>()},
                                                 {c2py::python_typename<double>()}},
                                                {c2py::python_typename<inchworm::frame_t>()});
//--------------------- module function table  -----------------------------

static PyMethodDef module_methods[] = {
   {"create_effective_hyb", (PyCFunction)c2py::pyfkw<_c2py_fun_0>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_0.c_str()},
   {"get_op_block_matrix", (PyCFunction)c2py::pyfkw<_c2py_fun_1>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_1.c_str()},
   {"make_ED_propagator", (PyCFunction)c2py::pyfkw<_c2py_fun_2>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_2.c_str()},
   {"make_bare_g_frame", (PyCFunction)c2py::pyfkw<_c2py_fun_3>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_3.c_str()},
   {"make_bare_u_frame", (PyCFunction)c2py::pyfkw<_c2py_fun_4>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_4.c_str()},
   {"make_g_frame_from_l_and_r", (PyCFunction)c2py::pyfkw<_c2py_fun_5>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_5.c_str()},
   {"partial_trace_bath", (PyCFunction)c2py::pyfkw<_c2py_fun_6>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_6.c_str()},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

//--------------------- module struct & init error definition ------------

//// module doc directly in the code or "" if not present...
/// Or mandatory ?
static struct PyModuleDef module_def = {PyModuleDef_HEAD_INIT,
                                        "ad_tools",                                /* name of module */
                                        R"RAWDOC(Inchworm atom_diag tools)RAWDOC", /* module documentation, may be NULL */
                                        -1, /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
                                        module_methods,
                                        NULL,
                                        NULL,
                                        NULL,
                                        NULL};

//--------------------- module init function -----------------------------

extern "C" __attribute__((visibility("default"))) PyObject *PyInit_ad_tools() {

  if (not c2py::check_python_version("ad_tools")) return NULL;

  // import numpy iff 'numpy/arrayobject.h' included
#ifdef Py_ARRAYOBJECT_H
  import_array();
#endif

  PyObject *m;

  if (PyType_Ready(&c2py::wrap_pytype<c2py::py_range>) < 0) return NULL;

  m = PyModule_Create(&module_def);
  if (m == NULL) return NULL;

  auto &conv_table = *c2py::conv_table_sptr.get();

  conv_table[std::type_index(typeid(c2py::py_range)).name()] = &c2py::wrap_pytype<c2py::py_range>;
#define _add_type(T, N) c2py::add_type_object_to_main<T>(N, m, conv_table)

#undef _add_type

  return m;
}
#endif
// CLAIR_WRAP_GEN
