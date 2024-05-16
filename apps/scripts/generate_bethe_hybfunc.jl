using ArgParse
using HDF5
using QuadGK

function spectral_density(omega, t)
    return 1 / (2 * t^2 * pi) * sqrt(4 * t^2 - omega^2)
end

function integrand(omega, tau, beta, t)
    return -1 * spectral_density(omega, t) / (exp(-beta * omega) + 1.0) * exp(-tau * omega)
end

function compute_delta(tau_grid, beta, t)
    delta_values = zeros(length(tau_grid))

    for itau in 1:length(tau_grid)
        tau = tau_grid[itau]
        delta_values[itau] = quadgk(omega->integrand(omega,tau,beta,t), -2*t, 2*t)[1]
    end

    return delta_values
end

function save_to_hdf5(file_path, bl_structure, tau_grid, data, spin_orb)
    h5open(file_path, "w") do file
        write(file, "bl_structure",  Int32.(bl_structure))
        write(file, "tau_grid", tau_grid)
        for ispin_orb in 0:(spin_orb-1)
            write(file, "data/$ispin_orb" * "_00", data)
        end    
end
end

function main(args)
    tau_grid = range(start=0., stop=args["beta"], length=Int(args["Nt"])) |> collect
    delta_values = compute_delta(tau_grid, args["beta"], args["t"])
    bl_structure = ones(Int,args["spin_orb"])
    save_to_hdf5(args["output_file"], bl_structure, tau_grid, delta_values, args["spin_orb"])
end

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "spin_orb"
            help = "number of spin-orbital"
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
        "output_file"
            help = "Output file path"
            arg_type = String
            default = "Delta_Bethe.h5"
    end
    return parse_args(s)
end

args = parse_commandline()
main(args)
