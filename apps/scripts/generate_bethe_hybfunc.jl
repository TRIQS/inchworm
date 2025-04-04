using ArgParse
using HDF5
using QuadGK
using Printf

function spectral_density(omega, t)
	return 1 / (2 * t^2 * pi) * sqrt(4 * t^2 - omega^2)
end

function integrand(omega, tau, beta, t)
	return -1 * spectral_density(omega, t) / (exp(-beta * omega) + 1.0) * exp(-tau * omega)
end

function compute_delta(tau_grid, beta, t, V)
	delta_values = zeros(length(tau_grid))

	for itau in 1:length(tau_grid)
		tau = tau_grid[itau]
		delta_values[itau] = quadgk(omega -> integrand(omega, tau, beta, t), -2 * t, 2 * t, atol = 1e-15, rtol = 1e-12)[1]
	end

	return V * V * delta_values
end

function save_to_hdf5(file_path, bl_structure, tau_grid, data, spin, orbital, off_diagonal_strength)
	h5open(file_path, "w") do file
		write(file, "bl_structure", Int32.(bl_structure))
		write(file, "tau_grid", tau_grid)
		for ispin in 1:spin
			for iorb_i in 1:orbital
				for iorb_j in 1:orbital
					if iorb_i == iorb_j
						prefactor = 1.0
					else
						prefactor = off_diagonal_strength
					end
					write(file, "data/$(ispin-1)_$(iorb_i-1)$(iorb_j-1)", prefactor * data)
				end
			end
		end
	end
end

function save_to_txt(txt_file_path, tau_grid, data, spin, orbital, off_diagonal_strength)
    spin_orb = spin * orbital
	open(txt_file_path, "w") do file
		for i in 1:length(tau_grid)
			for j in 1:spin_orb
                for k in 1:spin_orb
					i_0 = i-1
                    j_0 = j-1 
					k_0 = k-1
                    orbital_j = div(j_0, spin)
					spin_j = j_0%spin
					orbital_k = div(k_0, spin)
					spin_k = k_0%spin
					if spin_j != spin_k
						prefactor = 0.0
					else
						if orbital_j == orbital_k
							prefactor = 1.0
						else
							prefactor = off_diagonal_strength
						end
					end
					# Format the line with consistent spacing
					line = @sprintf(
						"%d %d %d %.16f %.1f\n",
						i_0, j_0, k_0, prefactor * data[i], 0
					)
					write(file, line)
				end
			end
		end
	end
end


function main(args)
	args["beta"] = args["beta"] * args["rescale"]
	V = 1.0
	V = V / args["rescale"]
	args["t"] = args["t"] / args["rescale"]
	tau_grid = range(start = 0.0, stop = args["beta"], length = Int(args["Nt"])) |> collect
	delta_values = compute_delta(tau_grid, args["beta"], args["t"], V)
	bl_structure = ones(Int, args["spin"]) * args["orbital"]
	save_to_hdf5(args["output_file"], bl_structure, tau_grid, delta_values, args["spin"], args["orbital"], args["off_diagonal_strength"])
	if args["save_txt"]
		txt_file_path = replace(args["output_file"], ".h5" => ".txt")
		save_to_txt(txt_file_path, tau_grid, delta_values, args["spin"], args["orbital"], args["off_diagonal_strength"])
	end
end

function parse_commandline()
	s = ArgParseSettings()
	@add_arg_table s begin
		"spin"
		help = "number of spin"
		arg_type = Int32
		required = true
		"orbital"
		help = "number of orbital"
		arg_type = Int32
		required = true
		"beta"
		help = "inverse temperature"
		arg_type = Float64
		required = true
		"Nt"
		help = "number of points for tau"
		arg_type = Float64
		required = true
		"t"
		help = "hopping"
		arg_type = Float64
		default = 1.0
		"off_diagonal_strength"
		help = "off-diagonal hybridization strength"
		arg_type = Float64
		default = 1.0
		"output_file"
		help = "Output file path"
		arg_type = String
		default = "Delta_Bethe.h5"
		"save_txt"
		help = "Flag to save data to a .txt file"
		arg_type = Bool
		default = false
		"rescale"
		help = "rescale beta to beta*rescale and V^2 to V^2/rescale^2"
		arg_type = Float64
		default = 1.0
	end
	return parse_args(s)
end

args = parse_commandline()
main(args)
