###############################################################################
#
# inchworm: A TRIQS based impurity solver
#
# Copyright (c) 2019--present, The Simons Foundation
#
# inchworm is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# inchworm is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# inchworm. If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################
from .solver_core import SolverCore, ConstrParamsT, SolveParamsT
from .version import is_delta_complex

from triqs.gfs import *
from triqs.gfs.tools import discretize_bath
from triqs.utility import mpi
from triqs.operators import *

import numpy

# === The SolverCore Wrapper

def mpi_print(*args, **kwargs):
    if mpi.is_master_node():
        print(*args, **kwargs)

class Solver(SolverCore):
    def __init__(self, beta, gf_struct, **kwargs):
        """
        Initialise the solver.

        Parameters
        ----------
        beta : scalar
               Inverse temperature.
        gf_struct : list of pairs [ [str,[int,...]], ...]
                    Structure of the Green's functions. It must be a
                    list of pairs, each containing the name of the
                    Green's function block as a string and a list of integer
                    indices.
                    For example: ``[ ['up', [0, 1, 2]], ['down', [0, 1, 2]] ]``.
        n_iw : integer [default=500]
               Number of Matsubara frequencies used for the Green's functions.
        n_tau : integer [default=5001]
               Number of imaginary time points used for the Green's functions.
        """
        # Initialise the core solver
        SolverCore.__init__(self, ConstrParamsT(beta=beta, gf_struct=gf_struct, **kwargs))

    def fit_bath(self, n_bath_sites_ED):
        if n_bath_sites_ED == 0:
            self.h_hyb_ED = Operator()
            self.Delta_ED = None
            return

        mpi_print("Finding discrete approximation for each hybridization block using {} sites ..".format(n_bath_sites_ED))
        V_mats, eps_lst, self.Delta_tau_ED = discretize_bath(self.Delta_tau, n_bath_sites_ED, eps0=3, V0=0.1, tol=1e-8, maxiter = 100000, cmplx=is_delta_complex, method='BFGS')

        gf_struct = self.constr_params["gf_struct"]
        bl_sizes = dict(gf_struct)
        bl_names = [bl_name for bl_name, bl_size in gf_struct]

        c_dag_vec = { b: numpy.matrix([[c_dag(b,o) for o in range(bl_sizes[b] + n_bath_sites_ED)]]) for b in bl_names }
        c_vec =     { b: numpy.matrix([[c(b,o)] for o in range(bl_sizes[b] + n_bath_sites_ED)]) for b in bl_names }

        self.h_hyb_ED = sum(c_dag_vec[b] * numpy.block([[numpy.zeros((bl_sizes[b], bl_sizes[b])), V_mats[n]], [V_mats[n].conj().T, numpy.diag(eps_lst[n])]]) * c_vec[b] for n, b in enumerate(bl_names))[0,0]
        mpi_print("", flush=True)

    def solve(self, **kwargs):
        """
        Solve the impurity problem.

        Parameters
        ----------
        params_kw : dict {'param':value} that is passed to the core solver.
                     The only two required parameters are
                        * `h_imp`: The local interaction Hamiltonian
                        * `n_cycles`: The number of Monte-Carlo cycles
                     For the other optional parameters see documentation.
                     Note that in this Python Wrapper the alpha-tensor is optional.
                     If not given, it will be constructed from the density matrix of
                     the SC Hartree Fock solution.
        """
        self.fit_bath(kwargs.get("n_bath_sites_ED",0))

        # Call the core solver's solve routine
        return SolverCore.solve(self, SolveParamsT(**kwargs))

    def solve_self_consistently(self, solve_params, *args):
        self.fit_bath(solve_params.get("n_bath_sites_ED",0))
        return SolverCore.solve_self_consistently(self, SolveParamsT(**solve_params), *args)

    def solve_inchworm(self, solve_params, *args):
        self.fit_bath(solve_params.get("n_bath_sites_ED",0))
        return SolverCore.solve_inchworm(self, SolveParamsT(**solve_params), *args)

    def solve_green(self, **kwargs):
        self.fit_bath(kwargs.get("n_bath_sites_ED",0))
        return SolverCore.solve_green(self, SolveParamsT(**kwargs))
