import argparse
import numpy as np
from scipy.integrate import quad
import h5py

def spectral_density(omega, t):
    return 1 / (2 * t**2 * np.pi) * np.sqrt(4 * t**2 - omega**2)

def integrand(omega, tau, beta, t):
    return -1 * spectral_density(omega, t) / (np.exp(-beta * omega) + 1.0) * np.exp(-tau * omega)

def compute_delta(tau_grid, beta, t):
    delta_values = np.zeros((2, 2, len(tau_grid)))
    
    for itau, tau in enumerate(tau_grid):
        delta_values[0, 0, itau] = quad(lambda w: integrand(w, tau, beta, t), -2 * t, 2 * t)[0]
        delta_values[1, 1, itau] = delta_values[0, 0, itau]
    
    return delta_values

def save_to_hdf5(file_path, bl_structure, tau_grid, data):
    with h5py.File(file_path, 'w') as f:
        f['bl_structure'] = np.array(bl_structure)
        f['tau_grid'] = tau_grid
        f['data/00'] = data[0, 0, :]
        f['data/11'] = data[1, 1, :]

def main(args):
    tau_grid = np.linspace(0, args.beta, int(args.Nt))

    delta_values = compute_delta(tau_grid, args.beta, args.t)

    save_to_hdf5(args.output_file, [1, 1], tau_grid, delta_values)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='compute and save hybridization function for the Bethe lattice to HDF5.')
    parser.add_argument('--beta', type=float, help='inverse temperature', required=True)
    parser.add_argument('--Nt', type=float, help='number of points for tau', required=True)
    parser.add_argument('--t', type=float, default=1.0, help='hopping')
    parser.add_argument('--output_file', type=str, default='Delta_Bethe.h5', help='Output file path')

    args = parser.parse_args()
    main(args)



