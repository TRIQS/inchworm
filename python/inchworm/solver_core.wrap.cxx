
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
#include <c2py/serialization/h5.hpp>

using c2py::operator""_a;

// ==================== enums =====================

// ==================== module classes =====================

// --------- class _c2py_cls_0 -----------
using _c2py_cls_0                                            = inchworm::fop_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_0>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_0> = "inchworm.solver_core.FopT";

static int synth_constructor_0(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(PyExc_RuntimeError, ("Error in constructing inchworm::fop_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_0> *)self)->_c = new _c2py_cls_0{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError, ("Error in constructing inchworm::fop_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_0> *)self)->_c);
  de("tau", self_c.tau, false);
  de("dag", self_c.dag, false);
  de("linear_index", self_c.linear_index, false);
  de("bl", self_c.bl, false);
  de("idx", self_c.idx, false);
  de("left_width", self_c.left_width, true);
  de("right_width", self_c.right_width, true);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_0> = synth_constructor_0;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_0> =
   c2py::replace_tags(R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
tau : {par_0}

dag : {par_1}

linear_index : {par_2}

bl : {par_3}

idx : {par_4}

left_width : {par_5}, default=100.0

right_width : {par_6}, default=100.0

)DOC",
                      "par",
                      {c2py::python_typename<double>(), c2py::python_typename<bool>(), c2py::python_typename<long>(), c2py::python_typename<long>(),
                       c2py::python_typename<long>(), c2py::python_typename<double>(), c2py::python_typename<double>()});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_0>[] = {

   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_0 = R"DOC(The imaginary time // FIXME Or time_pt ?)DOC";
constexpr auto _c2py_doc_member_1 = R"DOC(C (false) or Cdag (true) // FIXME Can we make this compile-time?)DOC";
constexpr auto _c2py_doc_member_2 = R"DOC(The linear index in the fundamental_operator_set)DOC";
constexpr auto _c2py_doc_member_3 = R"DOC(The block index)DOC";
constexpr auto _c2py_doc_member_4 = R"DOC(The orbital (or non-block) index)DOC";
constexpr auto _c2py_doc_member_5 = R"DOC(The left and right width of the insertion-move tau distribution)DOC";
constexpr auto _c2py_doc_member_6 = R"DOC()DOC";
static PyObject *prop_get_dict_0(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_0> *)self)->_c);
  c2py::pydict dic;
  dic["tau"]          = self_c.tau;
  dic["dag"]          = self_c.dag;
  dic["linear_index"] = self_c.linear_index;
  dic["bl"]           = self_c.bl;
  dic["idx"]          = self_c.idx;
  dic["left_width"]   = self_c.left_width;
  dic["right_width"]  = self_c.right_width;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_0>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_0::tau, _c2py_cls_0>("tau", _c2py_doc_member_0),
   c2py::getsetdef_from_member<&_c2py_cls_0::dag, _c2py_cls_0>("dag", _c2py_doc_member_1),
   c2py::getsetdef_from_member<&_c2py_cls_0::linear_index, _c2py_cls_0>("linear_index", _c2py_doc_member_2),
   c2py::getsetdef_from_member<&_c2py_cls_0::bl, _c2py_cls_0>("bl", _c2py_doc_member_3),
   c2py::getsetdef_from_member<&_c2py_cls_0::idx, _c2py_cls_0>("idx", _c2py_doc_member_4),
   c2py::getsetdef_from_member<&_c2py_cls_0::left_width, _c2py_cls_0>("left_width", _c2py_doc_member_5),
   c2py::getsetdef_from_member<&_c2py_cls_0::right_width, _c2py_cls_0>("right_width", _c2py_doc_member_6),
   {"__dict__", (getter)prop_get_dict_0, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_0> = R"DOC(Type representing a creation or annihilation operator at a fixed time)DOC"
   + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_0>;
// --------- class _c2py_cls_1 -----------
using _c2py_cls_1                                            = inchworm::qmc_results_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_1>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_1> = "inchworm.solver_core.QmcResultsT";
static const auto _c2py_init_0 = c2py::dispatcher_c_kw_t{c2py::c_constructor<_c2py_cls_1, const std::vector<int> &>("shape_of_frame")};
template <> constexpr initproc c2py::tp_init<_c2py_cls_1>    = c2py::pyfkw_constructor<_c2py_init_0>;
template <> const std::string c2py::tp_ctor_doc<_c2py_cls_1> = _c2py_init_0.doc(R"DOC()DOC");
// normalize
static auto const _c2py_fun_0 = c2py::dispatcher_f_kw_t{
   c2py::cmethod([](_c2py_cls_1 &self, inchworm::scalar_t normalization_cte) -> decltype(auto) { return self.normalize(normalization_cte); }, "self",
                 "normalization_cte")};

// print
static auto const _c2py_fun_1 = c2py::dispatcher_f_kw_t{
   c2py::cmethod([](_c2py_cls_1 &self, int verbosity) -> decltype(auto) { return self.print(verbosity); }, "self", "verbosity"_a = 4)};

static const auto _c2py_doc_0 = _c2py_fun_0.doc(R"DOC()DOC");
static const auto _c2py_doc_1 = _c2py_fun_1.doc(R"DOC()DOC");

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_1>[] = {
   {"normalize", (PyCFunction)c2py::pyfkw<_c2py_fun_0>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_0.c_str()},
   {"print", (PyCFunction)c2py::pyfkw<_c2py_fun_1>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_1.c_str()},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_7  = R"DOC()DOC";
constexpr auto _c2py_doc_member_8  = R"DOC()DOC";
constexpr auto _c2py_doc_member_9  = R"DOC()DOC";
constexpr auto _c2py_doc_member_10 = R"DOC()DOC";
constexpr auto _c2py_doc_member_11 = R"DOC()DOC";
constexpr auto _c2py_doc_member_12 = R"DOC()DOC";
constexpr auto _c2py_doc_member_13 = R"DOC()DOC";
constexpr auto _c2py_doc_member_14 = R"DOC()DOC";
constexpr auto _c2py_doc_member_15 = R"DOC()DOC";
constexpr auto _c2py_doc_member_16 = R"DOC()DOC";

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_1>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_1::frame, _c2py_cls_1>("frame", _c2py_doc_member_7),
   c2py::getsetdef_from_member<&_c2py_cls_1::weight_zeroth_order, _c2py_cls_1>("weight_zeroth_order", _c2py_doc_member_8),
   c2py::getsetdef_from_member<&_c2py_cls_1::errs_frame, _c2py_cls_1>("errs_frame", _c2py_doc_member_9),
   c2py::getsetdef_from_member<&_c2py_cls_1::frame_by_order, _c2py_cls_1>("frame_by_order", _c2py_doc_member_10),
   c2py::getsetdef_from_member<&_c2py_cls_1::err_frame_by_order, _c2py_cls_1>("err_frame_by_order", _c2py_doc_member_11),
   c2py::getsetdef_from_member<&_c2py_cls_1::order_histogram, _c2py_cls_1>("order_histogram", _c2py_doc_member_12),
   c2py::getsetdef_from_member<&_c2py_cls_1::average_sign, _c2py_cls_1>("average_sign", _c2py_doc_member_13),
   c2py::getsetdef_from_member<&_c2py_cls_1::average_order, _c2py_cls_1>("average_order", _c2py_doc_member_14),
   c2py::getsetdef_from_member<&_c2py_cls_1::status, _c2py_cls_1>("status", _c2py_doc_member_15),
   c2py::getsetdef_from_member<&_c2py_cls_1::auto_corr_time, _c2py_cls_1>("auto_corr_time", _c2py_doc_member_16),

   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <> const std::string c2py::tp_doc<_c2py_cls_1> = R"DOC()DOC" + c2py::tp_ctor_doc<_c2py_cls_1>;
// --------- class _c2py_cls_2 -----------
using _c2py_cls_2                                            = inchworm::constr_params_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_2>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_2> = "inchworm.solver_core.ConstrParamsT";

static int synth_constructor_1(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing inchworm::constr_params_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_2> *)self)->_c = new _c2py_cls_2{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError, ("Error in constructing inchworm::constr_params_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_2> *)self)->_c);
  de("n_tau", self_c.n_tau, true);
  de("n_tau_inch", self_c.n_tau_inch, true);
  de("n_tau_green", self_c.n_tau_green, true);
  de("n_iw", self_c.n_iw, true);
  de("beta", self_c.beta, false);
  de("gf_struct", self_c.gf_struct, false);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_2> = synth_constructor_1;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_2> =
   c2py::replace_tags(R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
beta : {par_0}

gf_struct : {par_1}

n_tau : {par_2}, default=10001

n_tau_inch : {par_3}, default=101

n_tau_green : {par_4}, default=101

n_iw : {par_5}, default=5

)DOC",
                      "par",
                      {c2py::python_typename<double>(), c2py::python_typename<triqs::gfs::gf_struct_t>(), c2py::python_typename<int>(),
                       c2py::python_typename<int>(), c2py::python_typename<int>(), c2py::python_typename<int>()});
// block_names
static auto const _c2py_fun_2 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_2 const &self) -> decltype(auto) { return self.block_names(); }, "self")};

// n_blocks
static auto const _c2py_fun_3 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_2 const &self) -> decltype(auto) { return self.n_blocks(); }, "self")};

static const auto _c2py_doc_2 = _c2py_fun_2.doc(R"DOC(
Names of block indeces for the Green function
)DOC");
static const auto _c2py_doc_3 = _c2py_fun_3.doc(R"DOC(
Number of block indeces for the Green function
)DOC");

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_2>[] = {
   {"block_names", (PyCFunction)c2py::pyfkw<_c2py_fun_2>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_2.c_str()},
   {"n_blocks", (PyCFunction)c2py::pyfkw<_c2py_fun_3>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_3.c_str()},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_17 = R"DOC(Number of tau points for the hybridization function)DOC";
constexpr auto _c2py_doc_member_18 = R"DOC(Number of tau points for the propagator)DOC";
constexpr auto _c2py_doc_member_19 = R"DOC(Number of tau points for the Green function)DOC";
constexpr auto _c2py_doc_member_20 = R"DOC(Number of Matsubara frequencies)DOC";
constexpr auto _c2py_doc_member_21 = R"DOC(Inverse temperature)DOC";
constexpr auto _c2py_doc_member_22 = R"DOC(Block structure of the gf)DOC";
static PyObject *prop_get_dict_1(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_2> *)self)->_c);
  c2py::pydict dic;
  dic["n_tau"]       = self_c.n_tau;
  dic["n_tau_inch"]  = self_c.n_tau_inch;
  dic["n_tau_green"] = self_c.n_tau_green;
  dic["n_iw"]        = self_c.n_iw;
  dic["beta"]        = self_c.beta;
  dic["gf_struct"]   = self_c.gf_struct;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_2>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_2::n_tau, _c2py_cls_2>("n_tau", _c2py_doc_member_17),
   c2py::getsetdef_from_member<&_c2py_cls_2::n_tau_inch, _c2py_cls_2>("n_tau_inch", _c2py_doc_member_18),
   c2py::getsetdef_from_member<&_c2py_cls_2::n_tau_green, _c2py_cls_2>("n_tau_green", _c2py_doc_member_19),
   c2py::getsetdef_from_member<&_c2py_cls_2::n_iw, _c2py_cls_2>("n_iw", _c2py_doc_member_20),
   c2py::getsetdef_from_member<&_c2py_cls_2::beta, _c2py_cls_2>("beta", _c2py_doc_member_21),
   c2py::getsetdef_from_member<&_c2py_cls_2::gf_struct, _c2py_cls_2>("gf_struct", _c2py_doc_member_22),
   {"__dict__", (getter)prop_get_dict_1, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_2> =
   R"DOC(The parameters for the solver construction)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_2>;
// --------- class _c2py_cls_3 -----------
using _c2py_cls_3                                            = inchworm::solve_params_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_3>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_3> = "inchworm.solver_core.SolveParamsT";

static int synth_constructor_2(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(PyExc_RuntimeError, ("Error in constructing inchworm::solve_params_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_3> *)self)->_c = new _c2py_cls_3{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError, ("Error in constructing inchworm::solve_params_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_3> *)self)->_c);
  de("h_imp", self_c.h_imp, false);
  de("partition_method", self_c.partition_method, true);
  de("quantum_numbers", self_c.quantum_numbers, true);
  de("n_bath_sites_ED", self_c.n_bath_sites_ED, true);
  de("max_prob_zeroth_order", self_c.max_prob_zeroth_order, true);
  de("n_tau_inch_stop", self_c.n_tau_inch_stop, true);
  de("n_cycles", self_c.n_cycles, false);
  de("length_cycle", self_c.length_cycle, true);
  de("n_warmup_cycles", self_c.n_warmup_cycles, true);
  de("n_callibration_cycles", self_c.n_callibration_cycles, true);
  de("max_order", self_c.max_order, true);
  de("random_seed", self_c.random_seed, true);
  de("random_name", self_c.random_name, true);
  de("use_double_insertion", self_c.use_double_insertion, true);
  de("max_time", self_c.max_time, true);
  de("verbosity", self_c.verbosity, true);
  de("measure_average_sign", self_c.measure_average_sign, true);
  de("measure_average_order", self_c.measure_average_order, true);
  de("measure_order_histogram", self_c.measure_order_histogram, true);
  de("measure_frame_by_order", self_c.measure_frame_by_order, true);
  de("post_process", self_c.post_process, true);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_3> = synth_constructor_2;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_3> = c2py::replace_tags(R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
h_imp : {par_0}

n_cycles : {par_1}

partition_method : {par_2}, default="automatic"

quantum_numbers : {par_3}, default={}

n_bath_sites_ED : {par_4}, default=0

max_prob_zeroth_order : {par_5}, default=0.5

n_tau_inch_stop : {par_6}, default={}

length_cycle : {par_7}, default={}

n_warmup_cycles : {par_8}, default=4096

n_callibration_cycles : {par_9}, default=4096

max_order : {par_10}, default={}

random_seed : {par_11}, default=34789 + 928374 * mpi::communicator().rank()

random_name : {par_12}, default=""

use_double_insertion : {par_13}, default=true

max_time : {par_14}, default=-1

verbosity : {par_15}, default== 0 ? 1 : 0

measure_average_sign : {par_16}, default=true

measure_average_order : {par_17}, default=true

measure_order_histogram : {par_18}, default=false

measure_frame_by_order : {par_19}, default=false

post_process : {par_20}, default=true

)DOC",
                                                                      "par",
                                                                      {c2py::python_typename<inchworm::many_body_op_t>(),
                                                                       c2py::python_typename<int>(),
                                                                       c2py::python_typename<std::string>(),
                                                                       c2py::python_typename<std::vector<inchworm::many_body_op_t>>(),
                                                                       c2py::python_typename<int>(),
                                                                       c2py::python_typename<double>(),
                                                                       c2py::python_typename<std::optional<int>>(),
                                                                       c2py::python_typename<std::optional<int>>(),
                                                                       c2py::python_typename<int>(),
                                                                       c2py::python_typename<int>(),
                                                                       c2py::python_typename<std::optional<int>>(),
                                                                       c2py::python_typename<int>(),
                                                                       c2py::python_typename<std::string>(),
                                                                       c2py::python_typename<int>(),
                                                                       c2py::python_typename<int>(),
                                                                       c2py::python_typename<int>(),
                                                                       c2py::python_typename<bool>(),
                                                                       c2py::python_typename<bool>(),
                                                                       c2py::python_typename<bool>(),
                                                                       c2py::python_typename<bool>(),
                                                                       c2py::python_typename<bool>()});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_3>[] = {

   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_23 = R"DOC(Impurity Hamiltonian)DOC";
constexpr auto _c2py_doc_member_24 = R"DOC(Partition method
type: str)DOC";
constexpr auto _c2py_doc_member_25 = R"DOC(Quantum numbers
type: list(Operator)
default: [Total Particle Number])DOC";
constexpr auto _c2py_doc_member_26 = R"DOC(Discrete number of sites per block in Delta_tau used to approximate the bath)DOC";
constexpr auto _c2py_doc_member_27 = R"DOC(Maximum probability of order zero sampling)DOC";
constexpr auto _c2py_doc_member_28 = R"DOC(Number of inchworm steps before interruption)DOC";
constexpr auto _c2py_doc_member_29 = R"DOC(Number of MC cycles)DOC";
constexpr auto _c2py_doc_member_30 = R"DOC(Length of a MC cycles, auto-deduce if not provided)DOC";
constexpr auto _c2py_doc_member_31 = R"DOC(Number of warmup cycles)DOC";
constexpr auto _c2py_doc_member_32 = R"DOC(Number of callibration cycles)DOC";
constexpr auto _c2py_doc_member_33 = R"DOC(The maximum order [optional])DOC";
constexpr auto _c2py_doc_member_34 = R"DOC(Random seed of the random generator)DOC";
constexpr auto _c2py_doc_member_35 = R"DOC(Name of the random generator)DOC";
constexpr auto _c2py_doc_member_36 = R"DOC(Use double insertion)DOC";
constexpr auto _c2py_doc_member_37 = R"DOC(Maximum running time in seconds (-1 : no limit))DOC";
constexpr auto _c2py_doc_member_38 = R"DOC(Verbosity)DOC";
constexpr auto _c2py_doc_member_39 = R"DOC(Measure the average sign)DOC";
constexpr auto _c2py_doc_member_40 = R"DOC(Measure the average perturbation order)DOC";
constexpr auto _c2py_doc_member_41 = R"DOC(Measure the perturbation order histogram)DOC";
constexpr auto _c2py_doc_member_42 = R"DOC(Measure the frame order-resolved)DOC";
constexpr auto _c2py_doc_member_43 = R"DOC(Perform post processing)DOC";
static PyObject *prop_get_dict_2(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_3> *)self)->_c);
  c2py::pydict dic;
  dic["h_imp"]                   = self_c.h_imp;
  dic["partition_method"]        = self_c.partition_method;
  dic["quantum_numbers"]         = self_c.quantum_numbers;
  dic["n_bath_sites_ED"]         = self_c.n_bath_sites_ED;
  dic["max_prob_zeroth_order"]   = self_c.max_prob_zeroth_order;
  dic["n_tau_inch_stop"]         = self_c.n_tau_inch_stop;
  dic["n_cycles"]                = self_c.n_cycles;
  dic["length_cycle"]            = self_c.length_cycle;
  dic["n_warmup_cycles"]         = self_c.n_warmup_cycles;
  dic["n_callibration_cycles"]   = self_c.n_callibration_cycles;
  dic["max_order"]               = self_c.max_order;
  dic["random_seed"]             = self_c.random_seed;
  dic["random_name"]             = self_c.random_name;
  dic["use_double_insertion"]    = self_c.use_double_insertion;
  dic["max_time"]                = self_c.max_time;
  dic["verbosity"]               = self_c.verbosity;
  dic["measure_average_sign"]    = self_c.measure_average_sign;
  dic["measure_average_order"]   = self_c.measure_average_order;
  dic["measure_order_histogram"] = self_c.measure_order_histogram;
  dic["measure_frame_by_order"]  = self_c.measure_frame_by_order;
  dic["post_process"]            = self_c.post_process;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_3>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_3::h_imp, _c2py_cls_3>("h_imp", _c2py_doc_member_23),
   c2py::getsetdef_from_member<&_c2py_cls_3::partition_method, _c2py_cls_3>("partition_method", _c2py_doc_member_24),
   c2py::getsetdef_from_member<&_c2py_cls_3::quantum_numbers, _c2py_cls_3>("quantum_numbers", _c2py_doc_member_25),
   c2py::getsetdef_from_member<&_c2py_cls_3::n_bath_sites_ED, _c2py_cls_3>("n_bath_sites_ED", _c2py_doc_member_26),
   c2py::getsetdef_from_member<&_c2py_cls_3::max_prob_zeroth_order, _c2py_cls_3>("max_prob_zeroth_order", _c2py_doc_member_27),
   c2py::getsetdef_from_member<&_c2py_cls_3::n_tau_inch_stop, _c2py_cls_3>("n_tau_inch_stop", _c2py_doc_member_28),
   c2py::getsetdef_from_member<&_c2py_cls_3::n_cycles, _c2py_cls_3>("n_cycles", _c2py_doc_member_29),
   c2py::getsetdef_from_member<&_c2py_cls_3::length_cycle, _c2py_cls_3>("length_cycle", _c2py_doc_member_30),
   c2py::getsetdef_from_member<&_c2py_cls_3::n_warmup_cycles, _c2py_cls_3>("n_warmup_cycles", _c2py_doc_member_31),
   c2py::getsetdef_from_member<&_c2py_cls_3::n_callibration_cycles, _c2py_cls_3>("n_callibration_cycles", _c2py_doc_member_32),
   c2py::getsetdef_from_member<&_c2py_cls_3::max_order, _c2py_cls_3>("max_order", _c2py_doc_member_33),
   c2py::getsetdef_from_member<&_c2py_cls_3::random_seed, _c2py_cls_3>("random_seed", _c2py_doc_member_34),
   c2py::getsetdef_from_member<&_c2py_cls_3::random_name, _c2py_cls_3>("random_name", _c2py_doc_member_35),
   c2py::getsetdef_from_member<&_c2py_cls_3::use_double_insertion, _c2py_cls_3>("use_double_insertion", _c2py_doc_member_36),
   c2py::getsetdef_from_member<&_c2py_cls_3::max_time, _c2py_cls_3>("max_time", _c2py_doc_member_37),
   c2py::getsetdef_from_member<&_c2py_cls_3::verbosity, _c2py_cls_3>("verbosity", _c2py_doc_member_38),
   c2py::getsetdef_from_member<&_c2py_cls_3::measure_average_sign, _c2py_cls_3>("measure_average_sign", _c2py_doc_member_39),
   c2py::getsetdef_from_member<&_c2py_cls_3::measure_average_order, _c2py_cls_3>("measure_average_order", _c2py_doc_member_40),
   c2py::getsetdef_from_member<&_c2py_cls_3::measure_order_histogram, _c2py_cls_3>("measure_order_histogram", _c2py_doc_member_41),
   c2py::getsetdef_from_member<&_c2py_cls_3::measure_frame_by_order, _c2py_cls_3>("measure_frame_by_order", _c2py_doc_member_42),
   c2py::getsetdef_from_member<&_c2py_cls_3::post_process, _c2py_cls_3>("post_process", _c2py_doc_member_43),
   {"__dict__", (getter)prop_get_dict_2, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_3> =
   R"DOC(The parameters for the solve function)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_3>;
// --------- class _c2py_cls_4 -----------
using _c2py_cls_4                                            = inchworm::solver_core;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_4>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_4> = "inchworm.solver_core.SolverCore";
static const auto _c2py_init_1 = c2py::dispatcher_c_kw_t{c2py::c_constructor<_c2py_cls_4, const inchworm::constr_params_t &>("constr_params_")};
template <> constexpr initproc c2py::tp_init<_c2py_cls_4> = c2py::pyfkw_constructor<_c2py_init_1>;
template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_4> = _c2py_init_1.doc(R"DOC(
Construct a INCHWORM solver

Parameters
----------
constr_params_ : {par_0}
   Set of parameters specific to the INCHWORM solver
)DOC",
                                                                    {{c2py::python_typename<const inchworm::constr_params_t &>()}});
// post_process
static auto const _c2py_fun_4 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_4 &self) -> decltype(auto) { return self.post_process(); }, "self")};

// solve
static auto const _c2py_fun_5 = c2py::dispatcher_f_kw_t{
   c2py::cmethod([](_c2py_cls_4 &self, const inchworm::solve_params_t &solve_params) -> decltype(auto) { return self.solve(solve_params); }, "self",
                 "solve_params")};

// solve_cthyb
static auto const _c2py_fun_6 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_4 &self, const inchworm::solve_params_t &solve_params,
                                            double tau_max) -> decltype(auto) { return self.solve_cthyb(solve_params, tau_max); },
                                         "self", "solve_params", "tau_max")};

// solve_green
static auto const _c2py_fun_7 = c2py::dispatcher_f_kw_t{
   c2py::cmethod([](_c2py_cls_4 &self, const inchworm::solve_params_t &solve_params) -> decltype(auto) { return self.solve_green(solve_params); },
                 "self", "solve_params")};

// solve_inchworm
static auto const _c2py_fun_8 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_4 &self, const inchworm::solve_params_t &solve_params,
                                            bool use_cthyb) -> decltype(auto) { return self.solve_inchworm(solve_params, use_cthyb); },
                                         "self", "solve_params", "use_cthyb"_a = false)};

// solve_self_consistently
static auto const _c2py_fun_9 = c2py::dispatcher_f_kw_t{
   c2py::cmethod([](_c2py_cls_4 &self, const inchworm::solve_params_t &solve_params, const inchworm::u_tau_t &u_tau_, double tau_split,
                    double tau_max) -> decltype(auto) { return self.solve_self_consistently(solve_params, u_tau_, tau_split, tau_max); },
                 "self", "solve_params", "u_tau_", "tau_split", "tau_max")};

static const auto _c2py_doc_4 = _c2py_fun_4.doc(R"DOC()DOC");
static const auto _c2py_doc_5 = _c2py_fun_5.doc(R"DOC(
Solve method that performs INCHWORM calculation

Parameters
----------
solve_params : {par_0}
   Set of parameters specific to the INCHWORM run
)DOC",
                                                {{c2py::python_typename<const inchworm::solve_params_t &>()}});
static const auto _c2py_doc_6 = _c2py_fun_6.doc(R"DOC()DOC");
static const auto _c2py_doc_7 = _c2py_fun_7.doc(R"DOC()DOC");
static const auto _c2py_doc_8 = _c2py_fun_8.doc(R"DOC()DOC");
static const auto _c2py_doc_9 = _c2py_fun_9.doc(R"DOC()DOC");

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_4>[] = {
   {"post_process", (PyCFunction)c2py::pyfkw<_c2py_fun_4>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_4.c_str()},
   {"solve", (PyCFunction)c2py::pyfkw<_c2py_fun_5>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_5.c_str()},
   {"solve_cthyb", (PyCFunction)c2py::pyfkw<_c2py_fun_6>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_6.c_str()},
   {"solve_green", (PyCFunction)c2py::pyfkw<_c2py_fun_7>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_7.c_str()},
   {"solve_inchworm", (PyCFunction)c2py::pyfkw<_c2py_fun_8>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_8.c_str()},
   {"solve_self_consistently", (PyCFunction)c2py::pyfkw<_c2py_fun_9>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_9.c_str()},
   {"__write_hdf5__", c2py::tpxx_write_h5<_c2py_cls_4>, METH_VARARGS, "  "},
   {"__getstate__", c2py::getstate_h5<_c2py_cls_4>, METH_NOARGS, ""},
   {"__setstate__", c2py::setstate_h5<_c2py_cls_4>, METH_O, ""},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_44 = R"DOC()DOC";
constexpr auto _c2py_doc_member_45 = R"DOC()DOC";
constexpr auto _c2py_doc_member_46 = R"DOC(Diagonalization of the local problem)DOC";
constexpr auto _c2py_doc_member_47 = R"DOC()DOC";
constexpr auto _c2py_doc_member_48 = R"DOC()DOC";
constexpr auto _c2py_doc_member_49 = R"DOC()DOC";
constexpr auto _c2py_doc_member_50 = R"DOC()DOC";
constexpr auto _c2py_doc_member_51 = R"DOC(The propagator in imaginary time)DOC";
constexpr auto _c2py_doc_member_52 = R"DOC(The zeroth order frame of the last sampling method)DOC";
constexpr auto _c2py_doc_member_53 = R"DOC()DOC";
constexpr auto _c2py_doc_member_54 = R"DOC(The vector of all creation operators)DOC";
constexpr auto _c2py_doc_member_55 = R"DOC(The vector of all annihilation operators)DOC";
constexpr auto _c2py_doc_member_56 = R"DOC(Greens function in imaginary time)DOC";
constexpr auto _c2py_doc_member_57 = R"DOC(Last obtained frame error block-resolved)DOC";
constexpr auto _c2py_doc_member_58 = R"DOC(Order-resolved propagator)DOC";
constexpr auto _c2py_doc_member_59 = R"DOC(Order-resolved propagator)DOC";
constexpr auto _c2py_doc_member_60 = R"DOC(Last obtained frame error order-resolved)DOC";
constexpr auto _c2py_doc_member_61 = R"DOC(Order histograms)DOC";

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_4>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_4::constr_params, _c2py_cls_4>("constr_params", _c2py_doc_member_44),
   c2py::getsetdef_from_member<&_c2py_cls_4::last_solve_params, _c2py_cls_4>("last_solve_params", _c2py_doc_member_45),
   c2py::getsetdef_from_member<&_c2py_cls_4::ad_imp, _c2py_cls_4>("ad_imp", _c2py_doc_member_46),
   c2py::getsetdef_from_member<&_c2py_cls_4::Delta_tau, _c2py_cls_4>("Delta_tau", _c2py_doc_member_47),
   c2py::getsetdef_from_member<&_c2py_cls_4::Delta_tau_ED, _c2py_cls_4>("Delta_tau_ED", _c2py_doc_member_48),
   c2py::getsetdef_from_member<&_c2py_cls_4::Delta_tau_tilde, _c2py_cls_4>("Delta_tau_tilde", _c2py_doc_member_49),
   c2py::getsetdef_from_member<&_c2py_cls_4::h_hyb_ED, _c2py_cls_4>("h_hyb_ED", _c2py_doc_member_50),
   c2py::getsetdef_from_member<&_c2py_cls_4::u_tau, _c2py_cls_4>("u_tau", _c2py_doc_member_51),
   c2py::getsetdef_from_member<&_c2py_cls_4::frame_zeroth_order, _c2py_cls_4>("frame_zeroth_order", _c2py_doc_member_52),
   c2py::getsetdef_from_member<&_c2py_cls_4::fops, _c2py_cls_4>("fops", _c2py_doc_member_53),
   c2py::getsetdef_from_member<&_c2py_cls_4::all_d_ops, _c2py_cls_4>("all_d_ops", _c2py_doc_member_54),
   c2py::getsetdef_from_member<&_c2py_cls_4::all_d_dag_ops, _c2py_cls_4>("all_d_dag_ops", _c2py_doc_member_55),
   c2py::getsetdef_from_member<&_c2py_cls_4::G_tau, _c2py_cls_4>("G_tau", _c2py_doc_member_56),
   c2py::getsetdef_from_member<&_c2py_cls_4::errs_frame, _c2py_cls_4>("errs_frame", _c2py_doc_member_57),
   c2py::getsetdef_from_member<&_c2py_cls_4::G_tau_by_order, _c2py_cls_4>("G_tau_by_order", _c2py_doc_member_58),
   c2py::getsetdef_from_member<&_c2py_cls_4::u_tau_by_order, _c2py_cls_4>("u_tau_by_order", _c2py_doc_member_59),
   c2py::getsetdef_from_member<&_c2py_cls_4::err_frame_by_order, _c2py_cls_4>("err_frame_by_order", _c2py_doc_member_60),
   c2py::getsetdef_from_member<&_c2py_cls_4::order_histograms, _c2py_cls_4>("order_histograms", _c2py_doc_member_61),

   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_4> = R"DOC(The Solver class)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_4>;

// ==================== module functions ====================

//--------------------- module function table  -----------------------------

static PyMethodDef module_methods[] = {
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

//--------------------- module struct & init error definition ------------

//// module doc directly in the code or "" if not present...
/// Or mandatory ?
static struct PyModuleDef module_def = {PyModuleDef_HEAD_INIT,
                                        "solver_core",                                   /* name of module */
                                        R"RAWDOC(The inchworm solve_core module)RAWDOC", /* module documentation, may be NULL */
                                        -1, /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
                                        module_methods,
                                        NULL,
                                        NULL,
                                        NULL,
                                        NULL};

//--------------------- module init function -----------------------------

extern "C" __attribute__((visibility("default"))) PyObject *PyInit_solver_core() {

  if (not c2py::check_python_version("solver_core")) return NULL;

  // import numpy iff 'numpy/arrayobject.h' included
#ifdef Py_ARRAYOBJECT_H
  import_array();
#endif

  PyObject *m;

  if (PyType_Ready(&c2py::wrap_pytype<c2py::py_range>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_0>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_1>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_2>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_3>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_4>) < 0) return NULL;

  m = PyModule_Create(&module_def);
  if (m == NULL) return NULL;

  auto &conv_table = *c2py::conv_table_sptr.get();

  conv_table[std::type_index(typeid(c2py::py_range)).name()] = &c2py::wrap_pytype<c2py::py_range>;
#define _add_type(T, N) c2py::add_type_object_to_main<T>(N, m, conv_table)
  _add_type(_c2py_cls_0, "FopT");
  _add_type(_c2py_cls_1, "QmcResultsT");
  _add_type(_c2py_cls_2, "ConstrParamsT");
  _add_type(_c2py_cls_3, "SolveParamsT");
  _add_type(_c2py_cls_4, "SolverCore");
#undef _add_type

  c2py::pyref module = c2py::pyref::module("h5.formats");
  if (not module) return nullptr;
  c2py::pyref register_class = module.attr("register_class");

  register_h5_type<_c2py_cls_4>(register_class);

  return m;
}
#endif
// CLAIR_WRAP_GEN
