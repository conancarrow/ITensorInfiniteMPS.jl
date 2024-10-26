using ITensors, ITensorMPS, ITensorInfiniteMPS, HDF5, Printf
include("../IFT_function.jl")

# Set parameter values
J = 1.0
h_values = collect(range(0, stop=2, length=201))
lambda_values = [0.01, 0.05, 0.1, 0.5]
bond_dimensions = [2, 3]

# Temporary storage for results to avoid recomputation
computed_results = Dict{Tuple{Int, Float64, Float64}, Dict}()

# Open the HDF5 file and write data
h5open("IFTData.h5", "w") do file
    for bond_dim in bond_dimensions
        for lambda in lambda_values
            for h in h_values
                try
                    # Compute and store result if not already in cache
                    if !haskey(computed_results, (bond_dim, lambda, h))
                        computed_results[(bond_dim, lambda, h)] = IFT_vumps(J, h, lambda, bond_dim)
                    end
                    result = computed_results[(bond_dim, lambda, h)]
                    
                    # Format paths with bond dimension, lambda, and h values
                    lambda_str = @sprintf("%.2f", lambda)
                    h_str = @sprintf("%.2f", h)
                    dataset_base_path = "/D$bond_dim/lambda_$lambda_str/h_$h_str"

                    # Write each key-value pair in the result dictionary to its own path
                    for (key, value) in result
                        dataset_path = joinpath(dataset_base_path, key)
                        write(file, dataset_path, value)
                    end
                catch e
                    println("Error encountered for bond_dim = $bond_dim, λ = $lambda, h = $h: $e")
                end
            end
        end
        
        # Write the "Vumps_Energy" vectors for each lambda across all h values for the current bond dimension
        for lambda in lambda_values
            try
                # Collect the Vumps_Energy values across all h values for the current lambda and bond_dim
                vumps_energy_vector = [computed_results[(bond_dim, lambda, h)]["Vumps_Energy"] for h in h_values]
                # Format lambda string with three digits
                lambda_str = @sprintf("%.2f", lambda)
                energy_path = "/D$bond_dim/Vumps_Energy/lambda_$lambda_str"
                write(file, energy_path, vumps_energy_vector)
            catch e
                println("Error encountered for bond_dim = $bond_dim, λ = $lambda while saving Vumps_Energy: $e")
            end
        end
    end
end

# Optionally, save the λ and h values as metadata for reference
h5write("IFTData.h5", "lambda_values", lambda_values)
h5write("IFTData.h5", "h_values", h_values)
h5write("IFTData.h5", "bond_dimensions", bond_dimensions)

