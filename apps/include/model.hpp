/*******************************************************************************
 *
 * inchworm: A TRIQS based impurity solver
 *
 * Copyright (c) 2019 The Simons foundation
 *   authors: Maxime Charlebois
 *
 * inchworm is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * inchworm is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * inchworm. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include <inchworm/atom_diag.hpp>
#include <inchworm/params.hpp>
#include <inchworm/u_frame.hpp>
#include <inchworm/util.hpp>
#include <inchworm/types.hpp>

#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <h5/h5.hpp>
#include <triqs/atom_diag/gf.hpp>

using namespace inchworm;
using mat_t = nda::matrix<double>;
using vec_t = nda::array<double, 1>;

// Prepare fundamental operator set
inline std::pair<fundamental_operator_set, std::vector<many_body_op_t>> make_fops(int n_site, int n_bath, int bath_offset, int n_spin) {
  fundamental_operator_set fops;
  std::vector<many_body_op_t> qn(1);
  for (int spin = 0; spin < n_spin; spin++)
    for (int i = 0; i < n_site; i++) {
      auto sp = ((spin == 0) ? "up" : "dn");
      fops.insert(sp, i);
      qn[0] += n(sp, i);
    }
  for (int spin = 0; spin < n_spin; spin++)
    for (int i = bath_offset; i < bath_offset + n_bath; i++) {
      auto sp = ((spin == 0) ? "up" : "dn");
      fops.insert(sp, i);
      qn[0] += n(sp, i);
    }
  return std::pair<fundamental_operator_set, std::vector<many_body_op_t>>(fops, qn);
}

// Exact non-interacting Green function of a single atomic level
inline double one_fermion(double tau, double eps, double beta) {
  if (eps >= 0) {
    return -std::exp(-tau * eps) / (1. + std::exp(-beta * eps));
  } else {
    return -std::exp((beta - tau) * eps) / (1. + std::exp(beta * eps));
  }
}

inline std::tuple<double, double, double, hyb_tau_t, atom_diag, u_tau_t, g_tau_t>
discrete_Hubbard_setup(int n_site, int n_bath, int n_spin, double U, double mu, double t, constr_params_t const &cp, mat_t const &theta,
                       vec_t const &eps, long n_tot, long n_tau_linear, int order_Chebyshev, std::vector<double> &grid_linear,
                       std::vector<double> &grid) {
  std::cout << "Discrete Hubbard model setup" << std::endl;

  // === Define fundamental operator sets

  // Full system
  auto [fops_tot, qn_tot] = make_fops(n_site, n_bath, n_site, n_spin);

  // Impurity
  auto [fops_imp, qn_imp] = make_fops(n_site, 0, n_site, n_spin);

  // Bath
  auto [fops_bath, qn_bath] = make_fops(0, n_bath, n_site, n_spin);

  // === Initialize Hamiltonians

  // h_imp: Hamiltonian of the impurity sites (n_site)
  many_body_operator h_imp, h_hyb, h_bath;
  for (int j = 0; j < n_site; j++) {
    h_imp -= mu * n("up", j);

    if (n_spin == 2) {
      h_imp -= mu * n("dn", j);
      h_imp += U * n("up", j) * n("dn", j);
    }
    for (int i = 0; i < n_site; i++) {
      if (i != j) {
        h_imp -= t * c_dag("up", i) * c("up", j);
        if (n_spin == 2) h_imp -= t * c_dag("dn", i) * c("dn", j);
      }
    }
  }

  //TEST: hidden fermion
  // double coeff = 1;
  // h_imp += coeff*c_dag("up", 0) * c("up", 2);
  // h_imp += coeff*c_dag("up", 2) * c("up", 0);
  // h_imp += coeff*c_dag("up", 1) * c("up", 2);
  // h_imp += coeff*c_dag("up", 2) * c("up", 1);

  // h_bath: Hamiltonian of the bath (n_site)
  for (int k = 0; k < n_bath; k++) {
    h_bath += eps(k) * n("up", k + n_site);
    if (n_spin == 2) h_bath += eps(k) * n("dn", k + n_site);
  }

  // h_hyb: Hamiltonian coupling the impurity and the bath
  for (int i = 0; i < n_site; i++) {
    for (int k = 0; k < n_bath; k++) {
      h_hyb += theta(i, k) * (c_dag("up", i) * c("up", k + n_site));
      h_hyb += theta(i, k) * (c_dag("up", k + n_site) * c("up", i));
      if (n_spin == 2) {
        h_hyb += theta(i, k) * (c_dag("dn", i) * c("dn", k + n_site));
        h_hyb += theta(i, k) * (c_dag("dn", k + n_site) * c("dn", i));
      }
    }
  }

  // === Define the 3 different atom_diag objects (ED calculation with Triqs)

  auto ad_tot  = inchworm::atom_diag(h_imp + h_bath + h_hyb, fops_tot);
  auto ad_imp  = inchworm::atom_diag(h_imp, fops_imp, create_effective_hyb(cp.gf_struct));
  auto ad_bath = inchworm::atom_diag(h_bath, fops_bath);

  // Calculate exact propagator
  u_tau_t u_tau = make_ED_propagator(ad_tot, ad_imp, ad_bath, cp.beta, n_tot, n_tau_linear, order_Chebyshev, grid_linear, grid);

  // Calculate exact Green function
  g_tau_t g_tau = real(atomic_g_tau(ad_tot, cp.beta, cp.gf_struct, cp.n_tau_green));

  auto Delta_tau = hyb_tau_t{{cp.beta, Fermion, cp.n_tau_hyb}, cp.gf_struct};

  // create hybridization:
  for (auto tau : Delta_tau[0].mesh()) {
    for (int block = 0; block < cp.gf_struct.size(); block++) {
      Delta_tau[block][tau] = 0.0;
      for (auto [i, j, n] : product_range(n_site, n_site, n_bath)) {
        Delta_tau[block][tau](i, j) += theta(i, n) * dagger(theta)(n, j) * one_fermion(tau, eps(n), cp.beta);
      }
    }
  }

  auto Z_bath              = partition_function(ad_bath, cp.beta);
  double Z_bath_correction = std::exp(-(ad_bath.get_gs_energy()) * cp.beta);
  double Z_imp_correction  = std::exp(-(ad_imp.get_gs_energy()) * cp.beta);

  // std::cout << "ad_bath.get_gs_energy() " << ad_bath.get_gs_energy() << std::endl;
  // std::cout << "ad_imp.get_gs_energy() " << ad_imp.get_gs_energy() << std::endl;
  //TEST: hidden fermion
  // for (auto tau : Delta_tau[0].mesh()) {
  //   for (int block = 0; block < cp.gf_struct.size(); block++) {
  //     for (auto [i, j] : product_range(2, 2)) {
  //       Delta_tau[block][tau](i, j) -= coeff*coeff*one_fermion(tau, mu, cp.beta);
  //     }
  //   }
  // }

  // std::cout << "Delta_tau[0](0,0) " << Delta_tau[0](0)(1, 0) << std::endl;
  // std::cout << "Delta_tau[0](0,0) " << Delta_tau[0](0)(0, 1) << std::endl;
  // //TEST: artificial off-diagonal hybridization
  // for (auto tau : Delta_tau[0].mesh()) {
  //   for (int block = 0; block < cp.gf_struct.size(); block++) {
  //     for (auto [i, j] : product_range(n_site, n_site)) {
  //       if (i > j){
  //       Delta_tau[block][tau](i, j)  += 100;
  //       }
  //       if(i < j){
  //       Delta_tau[block][tau](i, j)  -= 100;
  //       }
  //     }
  //   }
  // }
  // std::cout << "new Delta_tau[0](0,0) " << Delta_tau[0](0)(1, 0) << std::endl;
  // std::cout << "new Delta_tau[0](0,0) " << Delta_tau[0](0)(0, 1) << std::endl;
  return {Z_bath_correction, Z_imp_correction, Z_bath, Delta_tau, ad_imp, u_tau, g_tau};
}

inline atom_diag imp_Hubbard_setup(int n_site, int n_spin, double U, double mu, double t, constr_params_t const &cp) {
  std::cout << "imp_Hubbard_setup..." << std::endl;

  // === Define fundamental operator sets

  // Impurity
  auto [fops_imp, qn_imp] = make_fops(n_site, 0, n_site, n_spin);

  // === Initialize Hamiltonians

  // h_imp: Hamiltonian of the impurity sites (n_site)
  many_body_operator h_imp;
  for (int j = 0; j < n_site; j++) {
    h_imp -= mu * n("up", j);

    if (n_spin == 2) {
      h_imp -= mu * n("dn", j);
      h_imp += U * n("up", j) * n("dn", j);
    }
    for (int i = 0; i < n_site; i++) {
      if (i != j) {
        h_imp -= t * c_dag("up", i) * c("up", j);
        if (n_spin == 2) h_imp -= t * c_dag("dn", i) * c("dn", j);
      }
    }
  }

  // === Define the 1 different atom_diag objects (ED calculation with Triqs)

  auto ad_imp = inchworm::atom_diag(h_imp, fops_imp);

  // for (auto [bl, bl_size] : enumerate(ad_imp.get_subspace_dims())) {
  //   for (int i : range(bl_size)) {

  //     std::cout << "eigvalue:" << ad_imp.get_eigenvalue(bl, i) << std::endl; }
  // }
  // std::cout << " ground energy: " << ad_imp.get_gs_energy()<< std::endl;

  return ad_imp;
}

inline std::vector<std::vector<std::complex<double>>> read_hopping_file(int n_site, int n_spin, const std::string &file_name) {
  std::ifstream file(file_name);
  if (!file.is_open()) { throw std::runtime_error("Error[read_hopping_file]: Could not open file " + file_name); }
  std::vector<std::vector<std::complex<double>>> h; // i,j
  size_t n_flavors = n_site * n_spin;
  h.resize(n_flavors);
  int flavorI, flavorJ;
  double real, imag;
  for (size_t i = 0; i < n_flavors; ++i) {
    h[i].resize(n_flavors);
    for (size_t j = 0; j < n_flavors; ++j) {
      file >> flavorI >> flavorJ >> real >> imag;
      if (flavorI != i || flavorJ != j) { throw std::runtime_error("Error[read_hopping_file]: parsing hopping file: " + file_name); }
      h[i][j] = {real, imag};
    }
  }
  return h;
}

inline std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> read_interaction_file(int n_site, int n_spin,
                                                                                                      const std::string &file_name) {
  std::ifstream file(file_name);
  if (!file.is_open()) { throw std::runtime_error("Error[read_interaction_file]: Could not open file " + file_name); }

  size_t n_flavors = n_site * n_spin;
  // Initialize u to be a 4D vector of size n_flavors x n_flavors x n_flavors x
  // n_flavors with all elements set to 0
  std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> u;
  u.resize(n_flavors,
           std::vector<std::vector<std::vector<std::complex<double>>>>(
              n_flavors, std::vector<std::vector<std::complex<double>>>(n_flavors, std::vector<std::complex<double>>(n_flavors, {0, 0}))));

  size_t nNonzero = 0;
  file >> nNonzero;
  int flavorI, flavorJ, flavorK, flavorL;
  double real, imag;
  size_t count;
  for (size_t i = 0; i < nNonzero; ++i) {
    file >> count >> flavorI >> flavorJ >> flavorK >> flavorL >> real >> imag;
    std::cout << "read nonzero element [" << i << "]: count=" << count
              << ", flavorI=" << flavorI << ", flavorJ=" << flavorJ
              << ", flavorK=" << flavorK << ", flavorL=" << flavorL
              << " -> (" << real << "," << imag << ")" << std::endl;
    u[flavorI][flavorJ][flavorK][flavorL] = {real, imag};
    if (count != i) { throw std::runtime_error("Error[read_interaction_file]: parsing interaction file: " + file_name); }
    // if end of file is reached before nNonzero lines are read, throw an error
    if (file.eof() && i < nNonzero - 1) { throw std::runtime_error("Error[read_interaction_file]: parsing interaction file: " + file_name); }
  }
  return u;
}

inline atom_diag imp_step(int n_site, int n_spin, const std::string &interaction_file, double mu, const std::string &hopping_file,
                          constr_params_t const &cp) {
  std::cout << "imp_step..." << std::endl;

  auto interaction = read_interaction_file(n_site, n_spin, interaction_file);
  auto hopping     = read_hopping_file(n_site, n_spin, hopping_file);
  // === Define fundamental operator sets
  // Impurity
  auto [fops_imp, qn_imp] = make_fops(n_site, 0, n_site, n_spin);

  // === Initialize Hamiltonians
  //check if n_spin is 1 or 2 otherwise throw error
  if (n_spin != 1 && n_spin != 2) { throw std::runtime_error("Error[imp_step]: n_spin should be either 1 or 2"); }

  // h_imp: Hamiltonian of the impurity sites (n_site)
  many_body_operator h_imp;
  for (int j = 0; j < n_site; j++) {
    h_imp -= mu * n("up", j);
    if (n_spin == 2) { h_imp -= mu * n("dn", j); }
  }

  size_t n_site_spin = n_site * n_spin;
  // hopping part
  std::string spin_names[2] = {"up", "dn"};
  for (int i = 0; i < n_site_spin; i++) {
    for (int j = 0; j < n_site_spin; j++) {
      int spin_i = i % n_spin;
      int spin_j = j % n_spin;
      int site_i = i / n_spin;
      int site_j = j / n_spin;
      if (spin_i != spin_j) {
        if (hopping[i][j] != std::complex<double>(0, 0)) {
          throw std::runtime_error("Error[imp_step]: hopping matrix voilate spin conservation");
        } else {
          continue;
        }
      }
      h_imp += hopping[i][j].real() * c_dag(spin_names[spin_i], site_i) * c(spin_names[spin_j], site_j);
    }
  }

  //interaction part
  for (int i = 0; i < n_site_spin; i++) {
    for (int j = 0; j < n_site_spin; j++) {
      for (int k = 0; k < n_site_spin; k++) {
        for (int l = 0; l < n_site_spin; l++) {
          int spin_i = i % n_spin;
          int spin_j = j % n_spin;
          int spin_k = k % n_spin;
          int spin_l = l % n_spin;
          int site_i = i / n_spin;
          int site_j = j / n_spin;
          int site_k = k / n_spin;
          int site_l = l / n_spin;

          if (spin_i + spin_j != spin_k + spin_l) {
            if (interaction[i][j][k][l] != std::complex<double>(0, 0)) {
              throw std::runtime_error("Error[imp_step]: interaction matrix voilate spin conservation");
            } else {
              continue;
            }
          }
          h_imp += interaction[i][j][k][l].real() * c_dag(spin_names[spin_i], site_i) * c_dag(spin_names[spin_j], site_j)
             * c(spin_names[spin_k], site_k) * c(spin_names[spin_l], site_l);
        }
      }
    }
  }

  // === Define the 1 different atom_diag objects (ED calculation with Triqs)

  auto ad_imp = inchworm::atom_diag(h_imp, fops_imp);

  return ad_imp;
}


inline std::tuple<double, double, double, hyb_tau_t, atom_diag, u_tau_t, g_tau_t>
discrete_setup(int n_site, int n_bath, int n_spin, const std::string &interaction_file, double mu, const std::string &hopping_file, constr_params_t const &cp, mat_t const &theta,
                       vec_t const &eps, long n_tot, long n_tau_linear, int order_Chebyshev, std::vector<double> &grid_linear,
                       std::vector<double> &grid) {
  std::cout << "discrete_setup..." << std::endl;
  // === Define fundamental operator sets

  // Full system
  auto [fops_tot, qn_tot] = make_fops(n_site, n_bath, n_site, n_spin);

  // Impurity
  auto [fops_imp, qn_imp] = make_fops(n_site, 0, n_site, n_spin);

  // Bath
  auto [fops_bath, qn_bath] = make_fops(0, n_bath, n_site, n_spin);

  // === Initialize Hamiltonians

  // h_imp: Hamiltonian of the impurity sites (n_site)
  many_body_operator h_imp, h_hyb, h_bath;
  auto interaction = read_interaction_file(n_site, n_spin, interaction_file);
  auto hopping     = read_hopping_file(n_site, n_spin, hopping_file);
  if (n_spin != 1 && n_spin != 2) { throw std::runtime_error("Error[imp_step]: n_spin should be either 1 or 2"); }
  for (int j = 0; j < n_site; j++) {
    h_imp -= mu * n("up", j);
    if (n_spin == 2) { h_imp -= mu * n("dn", j); }
  }

  size_t n_site_spin = n_site * n_spin;
  // hopping part
  std::string spin_names[2] = {"up", "dn"};
  for (int i = 0; i < n_site_spin; i++) {
    for (int j = 0; j < n_site_spin; j++) {
      int spin_i = i % n_spin;
      int spin_j = j % n_spin;
      int site_i = i / n_spin;
      int site_j = j / n_spin;
      if (spin_i != spin_j) {
        if (hopping[i][j] != std::complex<double>(0, 0)) {
          throw std::runtime_error("Error[imp_step]: hopping matrix voilate spin conservation");
        } else {
          continue;
        }
      }
      h_imp += hopping[i][j].real() * c_dag(spin_names[spin_i], site_i) * c(spin_names[spin_j], site_j);
    }
  }

  //interaction part
  for (int i = 0; i < n_site_spin; i++) {
    for (int j = 0; j < n_site_spin; j++) {
      for (int k = 0; k < n_site_spin; k++) {
        for (int l = 0; l < n_site_spin; l++) {
          int spin_i = i % n_spin;
          int spin_j = j % n_spin;
          int spin_k = k % n_spin;
          int spin_l = l % n_spin;
          int site_i = i / n_spin;
          int site_j = j / n_spin;
          int site_k = k / n_spin;
          int site_l = l / n_spin;

          if (spin_i + spin_j != spin_k + spin_l) {
            if (interaction[i][j][k][l] != std::complex<double>(0, 0)) {
              throw std::runtime_error("Error[imp_step]: interaction matrix voilate spin conservation");
            } else {
              continue;
            }
          }
          h_imp += interaction[i][j][k][l].real() / 2.0 * c_dag(spin_names[spin_i], site_i) * c_dag(spin_names[spin_j], site_j)
             * c(spin_names[spin_k], site_k) * c(spin_names[spin_l], site_l);
        }
      }
    }
  }

  // h_bath: Hamiltonian of the bath (n_site)
  for (int k = 0; k < n_bath; k++) {
    h_bath += eps(k) * n("up", k + n_site);
    if (n_spin == 2) h_bath += eps(k) * n("dn", k + n_site);
  }

  // h_hyb: Hamiltonian coupling the impurity and the bath
  for (int i = 0; i < n_site; i++) {
    for (int k = 0; k < n_bath; k++) {
      h_hyb += theta(i, k) * (c_dag("up", i) * c("up", k + n_site));
      h_hyb += theta(i, k) * (c_dag("up", k + n_site) * c("up", i));
      if (n_spin == 2) {
        h_hyb += theta(i, k) * (c_dag("dn", i) * c("dn", k + n_site));
        h_hyb += theta(i, k) * (c_dag("dn", k + n_site) * c("dn", i));
      }
    }
  }

  // === Define the 3 different atom_diag objects (ED calculation with Triqs)

  auto ad_tot  = inchworm::atom_diag(h_imp + h_bath + h_hyb, fops_tot);
  auto ad_imp  = inchworm::atom_diag(h_imp, fops_imp, create_effective_hyb(cp.gf_struct));
  auto ad_bath = inchworm::atom_diag(h_bath, fops_bath);

  // Calculate exact propagator
  u_tau_t u_tau = make_ED_propagator(ad_tot, ad_imp, ad_bath, cp.beta, n_tot, n_tau_linear, order_Chebyshev, grid_linear, grid);

  // Calculate exact Green function
  g_tau_t g_tau = real(atomic_g_tau(ad_tot, cp.beta, cp.gf_struct, cp.n_tau_green));

  auto Delta_tau = hyb_tau_t{{cp.beta, Fermion, cp.n_tau_hyb}, cp.gf_struct};

  // create hybridization:
  for (auto tau : Delta_tau[0].mesh()) {
    for (int block = 0; block < cp.gf_struct.size(); block++) {
      Delta_tau[block][tau] = 0.0;
      for (auto [i, j, n] : product_range(n_site, n_site, n_bath)) {
        Delta_tau[block][tau](i, j) += theta(i, n) * dagger(theta)(n, j) * one_fermion(tau, eps(n), cp.beta);
      }
    }
  }

  auto Z_bath              = partition_function(ad_bath, cp.beta);
  double Z_bath_correction = std::exp(-(ad_bath.get_gs_energy()) * cp.beta);
  double Z_imp_correction  = std::exp(-(ad_imp.get_gs_energy()) * cp.beta);

  return {Z_bath_correction, Z_imp_correction, Z_bath, Delta_tau, ad_imp, u_tau, g_tau};
}
