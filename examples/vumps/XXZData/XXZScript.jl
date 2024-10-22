using ITensors, ITensorMPS, ITensorInfiniteMPS, HDF5

# Define the range of Delta (100 equally-spaced values from 0 to 4)
deltas = range(0, stop=4, length=201)

include("../xxz_function.jl")

# Iterate over each Delta value, compute xxz_vumps(x), and save to file
for Delta in deltas
    result = xxz_vumps(Delta)
    
    # Extract Delta and AL_array from the dictionary
    array_data = result["AL_Array"]
    delta_value = result["Delta"]
    
    # Create the filename using the Delta value
    filename = "XXZ_$(delta_value)_AL_Array.h5"
    
    # Open file and write the array to the file
    h5write(filename, "AL_Array", array_data)
    h5write(filename, "Delta", delta_value)
end

println("All files written successfully!")

