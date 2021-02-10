# Generated automatically using the command :
# c++2py ../../c++/inchworm/solver_core.hpp -p --members_read_only -N inchworm -a inchworm -m solver_core -o solver_core --moduledoc="The inchworm solve_core module" -C triqs -C nda_py --cxxflags="-std=c++20" --only="solver_core qmc_results_t"
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "solver_core", doc = r"The inchworm solve_core module", app_name = "inchworm")

# Imports
module.add_imports(*['triqs.atom_diag', 'triqs.gf', 'triqs.operators', 'h5._h5py'])

# Add here all includes
module.add_include("inchworm/solver_core.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/optional.hpp>
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


# The class qmc_results_t
c = class_(
        py_type = "QmcResultsT",  # name of the python class
        c_type = "inchworm::qmc_results_t",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_member(c_name = "frame",
             c_type = "inchworm::frame_t",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "frame_0th_order",
             c_type = "inchworm::frame_t",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "frame_by_order",
             c_type = "std::vector<frame_t>",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "order_histogram",
             c_type = "std::vector<double>",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "average_order",
             c_type = "double",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "status",
             c_type = "int",
             read_only= True,
             doc = r"""""")

c.add_constructor("""(std::vector<int> shape_of_frame)""", doc = r"""""")

c.add_method("""void normalize (double normalization_cte)""",
             doc = r"""""")

c.add_method("""void print (int verbosity = 4)""",
             doc = r"""""")

module.add_class(c)

# The class solver_core
c = class_(
        py_type = "SolverCore",  # name of the python class
        c_type = "inchworm::solver_core",   # name of the C++ class
        doc = r"""The Solver class""",   # doc of the C++ class
        hdf5 = True,
)

c.add_member(c_name = "G_tau",
             c_type = "inchworm::g_tau_t",
             read_only= True,
             doc = r"""Greens function in imaginary time""")

c.add_member(c_name = "order_histograms",
             c_type = "std::vector<std::vector<double>>",
             read_only= True,
             doc = r"""Order histograms""")

c.add_member(c_name = "constr_params",
             c_type = "inchworm::constr_params_t",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "last_solve_params",
             c_type = "std::optional<solve_params_t>",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "ad_imp",
             c_type = "inchworm::atom_diag",
             read_only= True,
             doc = r"""Diagonalization of the local problem""")

c.add_member(c_name = "Delta_tau",
             c_type = "inchworm::h_tau_t",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "u_tau",
             c_type = "inchworm::u_tau_t",
             read_only= True,
             doc = r"""The propagator in imaginary time""")

c.add_constructor("""(**inchworm::constr_params_t)""", doc = r"""Construct a INCHWORM solver



+----------------+-----------------------------------+---------+------------------------------------------------------+
| Parameter Name | Type                              | Default | Documentation                                        |
+================+===================================+=========+======================================================+
| n_tau          | int                               | 101     | Number of tau points for the hybridization function  |
+----------------+-----------------------------------+---------+------------------------------------------------------+
| n_tau_inch     | int                               | 101     | Number of tau points for the propagator              |
+----------------+-----------------------------------+---------+------------------------------------------------------+
| n_tau_green    | int                               | 101     | Number of tau points for the Green function          |
+----------------+-----------------------------------+---------+------------------------------------------------------+
| n_iw           | int                               | 5       | Number of Matsubara frequencies                      |
+----------------+-----------------------------------+---------+------------------------------------------------------+
| beta           | double                            | --      | Inverse temperature                                  |
+----------------+-----------------------------------+---------+------------------------------------------------------+
| gf_struct      | triqs::hilbert_space::gf_struct_t | --      | Block structure of the gf                            |
+----------------+-----------------------------------+---------+------------------------------------------------------+
""")

c.add_method("""void solve (**inchworm::solve_params_t)""",
             doc = r"""Solve method that performs INCHWORM calculation



+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| Parameter Name          | Type                                 | Default                                 | Documentation                                       |
+=========================+======================================+=========================================+=====================================================+
| h_imp                   | triqs::operators::many_body_operator | --                                      | Impurity Hamiltonian                                |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| partition_method        | str                                  | "automatic"                             | Partition method                                    |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| quantum_numbers         | list(Operator)                       | [Total Particle Number]                 | Quantum numbers                                     |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| max_prob_zeroth_order   | double                               | 1.0                                     | Maximum probability of order zero sampling          |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| n_cycles                | int                                  | --                                      | Number of MC cycles                                 |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| length_cycle            | std::optional<int>                   | {}                                      | Length of a MC cycles, auto-deduce if not provided  |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| n_warmup_cycles         | int                                  | 1000                                    | Number of warmup cycles                             |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| n_callibration_cycles   | int                                  | 1000                                    | Number of callibration cycles                       |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| max_order               | std::optional<int>                   | {}                                      | The maximum order [optional]                        |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| random_seed             | int                                  | 34789+928374*mpi::communicator().rank() | Random seed of the random generator                 |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| random_name             | std::string                          | ""                                      | Name of the random generator                        |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| use_double_insertion    | int                                  | true                                    | Use double insertion                                |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| max_time                | int                                  | -1                                      | Maximum running time in seconds (-1 : no limit)     |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| verbosity               | int                                  | mpi::communicator().rank()==0?1:0       | Verbosity                                           |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| measure_average_order   | bool                                 | true                                    | Measure the average perturbation order              |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| measure_order_histogram | bool                                 | false                                   | Measure the average perturbation order              |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| measure_frame_by_order  | bool                                 | false                                   | Measure the average perturbation order              |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| post_process            | bool                                 | true                                    | Perform post processing                             |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
""")

c.add_method("""inchworm::qmc_results_t solve_cthyb (inchworm::solve_params_t solve_params, double tau_max)""",
             doc = r"""""")

c.add_method("""inchworm::qmc_results_t solve_self_consistently (inchworm::solve_params_t solve_params, inchworm::u_tau_t u_tau_, double tau_split, double tau_max)""",
             doc = r"""""")

c.add_method("""void solve_inchworm (inchworm::solve_params_t solve_params, bool use_cthyb = false)""",
             doc = r"""""")

c.add_method("""void solve_green (**inchworm::solve_params_t)""",
             doc = r"""



+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| Parameter Name          | Type                                 | Default                                 | Documentation                                       |
+=========================+======================================+=========================================+=====================================================+
| h_imp                   | triqs::operators::many_body_operator | --                                      | Impurity Hamiltonian                                |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| partition_method        | str                                  | "automatic"                             | Partition method                                    |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| quantum_numbers         | list(Operator)                       | [Total Particle Number]                 | Quantum numbers                                     |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| max_prob_zeroth_order   | double                               | 1.0                                     | Maximum probability of order zero sampling          |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| n_cycles                | int                                  | --                                      | Number of MC cycles                                 |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| length_cycle            | std::optional<int>                   | {}                                      | Length of a MC cycles, auto-deduce if not provided  |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| n_warmup_cycles         | int                                  | 1000                                    | Number of warmup cycles                             |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| n_callibration_cycles   | int                                  | 1000                                    | Number of callibration cycles                       |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| max_order               | std::optional<int>                   | {}                                      | The maximum order [optional]                        |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| random_seed             | int                                  | 34789+928374*mpi::communicator().rank() | Random seed of the random generator                 |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| random_name             | std::string                          | ""                                      | Name of the random generator                        |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| use_double_insertion    | int                                  | true                                    | Use double insertion                                |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| max_time                | int                                  | -1                                      | Maximum running time in seconds (-1 : no limit)     |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| verbosity               | int                                  | mpi::communicator().rank()==0?1:0       | Verbosity                                           |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| measure_average_order   | bool                                 | true                                    | Measure the average perturbation order              |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| measure_order_histogram | bool                                 | false                                   | Measure the average perturbation order              |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| measure_frame_by_order  | bool                                 | false                                   | Measure the average perturbation order              |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
| post_process            | bool                                 | true                                    | Perform post processing                             |
+-------------------------+--------------------------------------+-----------------------------------------+-----------------------------------------------------+
""")

c.add_method("""std::string hdf5_format ()""",
             is_static = True,
             doc = r"""""")

c.add_property(name = "post_process",
               getter = cfunction("void post_process ()"),
               doc = r"""""")

module.add_class(c)


# Converter for solve_params_t
c = converter_(
        c_type = "inchworm::solve_params_t",
        doc = r"""The parameters for the solve function""",
)
c.add_member(c_name = "h_imp",
             c_type = "triqs::operators::many_body_operator",
             initializer = """  """,
             doc = r"""Impurity Hamiltonian""")

c.add_member(c_name = "partition_method",
             c_type = "std::string",
             initializer = """ "automatic" """,
             doc = r"""Partition method
     type: str""")

c.add_member(c_name = "quantum_numbers",
             c_type = "std::vector<many_body_op_t>",
             initializer = """ {} """,
             doc = r"""Quantum numbers
     type: list(Operator)
     default: [Total Particle Number]""")

c.add_member(c_name = "max_prob_zeroth_order",
             c_type = "double",
             initializer = """ 1.0 """,
             doc = r"""Maximum probability of order zero sampling""")

c.add_member(c_name = "n_cycles",
             c_type = "int",
             initializer = """  """,
             doc = r"""Number of MC cycles""")

c.add_member(c_name = "length_cycle",
             c_type = "std::optional<int>",
             initializer = """ {} """,
             doc = r"""Length of a MC cycles, auto-deduce if not provided""")

c.add_member(c_name = "n_warmup_cycles",
             c_type = "int",
             initializer = """ 1000 """,
             doc = r"""Number of warmup cycles""")

c.add_member(c_name = "n_callibration_cycles",
             c_type = "int",
             initializer = """ 1000 """,
             doc = r"""Number of callibration cycles""")

c.add_member(c_name = "max_order",
             c_type = "std::optional<int>",
             initializer = """ {} """,
             doc = r"""The maximum order [optional]""")

c.add_member(c_name = "random_seed",
             c_type = "int",
             initializer = """ 34789+928374*mpi::communicator().rank() """,
             doc = r"""Random seed of the random generator""")

c.add_member(c_name = "random_name",
             c_type = "std::string",
             initializer = """ "" """,
             doc = r"""Name of the random generator""")

c.add_member(c_name = "use_double_insertion",
             c_type = "int",
             initializer = """ true """,
             doc = r"""Use double insertion""")

c.add_member(c_name = "max_time",
             c_type = "int",
             initializer = """ -1 """,
             doc = r"""Maximum running time in seconds (-1 : no limit)""")

c.add_member(c_name = "verbosity",
             c_type = "int",
             initializer = """ mpi::communicator().rank()==0?1:0 """,
             doc = r"""Verbosity""")

c.add_member(c_name = "measure_average_order",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Measure the average perturbation order""")

c.add_member(c_name = "measure_order_histogram",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure the average perturbation order""")

c.add_member(c_name = "measure_frame_by_order",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure the average perturbation order""")

c.add_member(c_name = "post_process",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Perform post processing""")

module.add_converter(c)


# Converter for constr_params_t
c = converter_(
        c_type = "inchworm::constr_params_t",
        doc = r"""The parameters for the solver construction""",
)
c.add_member(c_name = "n_tau",
             c_type = "int",
             initializer = """ 101 """,
             doc = r"""Number of tau points for the hybridization function""")

c.add_member(c_name = "n_tau_inch",
             c_type = "int",
             initializer = """ 101 """,
             doc = r"""Number of tau points for the propagator""")

c.add_member(c_name = "n_tau_green",
             c_type = "int",
             initializer = """ 101 """,
             doc = r"""Number of tau points for the Green function""")

c.add_member(c_name = "n_iw",
             c_type = "int",
             initializer = """ 5 """,
             doc = r"""Number of Matsubara frequencies""")

c.add_member(c_name = "beta",
             c_type = "double",
             initializer = """  """,
             doc = r"""Inverse temperature""")

c.add_member(c_name = "gf_struct",
             c_type = "triqs::hilbert_space::gf_struct_t",
             initializer = """  """,
             doc = r"""Block structure of the gf""")

module.add_converter(c)


module.generate_code()
