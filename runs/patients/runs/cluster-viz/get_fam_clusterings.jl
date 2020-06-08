include("../../../sim-study/configs/test-sim-6-8-2/tsne/tsne_imports.jl")

function get_output_dir(path_to_output)
  joinpath(splitdir(path_to_output)[1], "output.bson")
end

### Main ###

# Results directory
results_dir = "$(ENV["SCRATCH_DIR"])/cytof/results/repfam/patients-data/run1"

# Where to put results
csvdir = "img"

# Loop through different settings
for phi in (0, 1, 10, 25)
  # Get simulation name
  simname = "phi$(phi)"
  println("Processing: ", simname)

  # Get path to simuldate data.
  path_to_output = joinpath(results_dir, simname, "output.bson")

  # Get output 
  output = BSON.load(path_to_output)
  
  # Get relevant parameters
  samples =  [s[1] for s in output[:samples][1]]

  # Get best lam (for each sample)
  best_lam = vcat(get_best_lam(samples)...)

  # Write output to file.
  open("$(csvdir)/fam-$(simname)-clusterings.csv", "w") do io
    writedlm(io, best_lam)
  end
end
