using ITensors, ITensorMPS, ITensorInfiniteMPS, HDF5
include("../IFT_function.jl")

# Set parameter values
J = 1.0
h_values = collect(range(0, stop=2, length=201))
lambda_values = [0.01, 0.05, 0.1, 0.5]

# Temporary storage for results to avoid recomputation
computed_results = Dict{Tuple{Float64, Float64}, Dict}()

# Open the HDF5 file and write data
h5open("D3_Data.h5", "w") do file
    for lambda in lambda_values
        for h in h_values
            try
                # Compute and store result if not already in cache
                if !haskey(computed_results, (lambda, h))
                    computed_results[(lambda, h)] = IFT_vumps(J, h, lambda)
                end
                result = computed_results[(lambda, h)]

                # Write each key-value pair in the result dictionary to its own path
                for (key, value) in result
                    dataset_path = "/lambda_$lambda/h_$h/$key"
                    write(file, dataset_path, value)
                end
            catch e
                println("Error encountered for λ = $lambda, h = $h: $e")
            end
        end
    end

    # Write the "Vumps_Energy" vectors for each lambda across all h values
    for lambda in lambda_values
        try
            # Collect the Vumps_Energy values across all h values for the current lambda
            vumps_energy_vector = [computed_results[(lambda, h)]["Vumps_Energy"] for h in h_values]
            # Save this vector as a single dataset
            energy_path = "/Vumps_Energy/lambda_$lambda"
            write(file, energy_path, vumps_energy_vector)
        catch e
            println("Error encountered for λ = $lambda while saving Vumps_Energy: $e")
        end
    end
end

# Optionally, save the λ and h values as metadata for reference
h5write("D3_Data.h5", "lambda_values", lambda_values)
h5write("D3_Data.h5", "h_values", h_values)

