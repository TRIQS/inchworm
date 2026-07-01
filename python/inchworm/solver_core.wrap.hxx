#include <c2py/c2py.hpp>

#ifndef C2PY_HXX_DECLARATION_solver_core_GUARDS
#define C2PY_HXX_DECLARATION_solver_core_GUARDS
template <> constexpr bool c2py::is_wrapped<inchworm::fop_t>               = true;
template <> inline constexpr auto c2py::tp_name<inchworm::fop_t>           = "inchworm.solver_core.FopT";
template <> constexpr bool c2py::is_wrapped<inchworm::qmc_results_t>       = true;
template <> inline constexpr auto c2py::tp_name<inchworm::qmc_results_t>   = "inchworm.solver_core.QmcResultsT";
template <> constexpr bool c2py::is_wrapped<inchworm::constr_params_t>     = true;
template <> inline constexpr auto c2py::tp_name<inchworm::constr_params_t> = "inchworm.solver_core.ConstrParamsT";
template <> constexpr bool c2py::is_wrapped<inchworm::solve_params_t>      = true;
template <> inline constexpr auto c2py::tp_name<inchworm::solve_params_t>  = "inchworm.solver_core.SolveParamsT";
template <> constexpr bool c2py::is_wrapped<inchworm::solver_core>         = true;
template <> inline constexpr auto c2py::tp_name<inchworm::solver_core>     = "inchworm.solver_core.SolverCore";
#endif