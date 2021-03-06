include("tsne_imports.jl")

get_output_path(path_to_output) = let
  joinpath(splitdir(path_to_output)[1], "output.bson")
end

results_dir = "$(ENV["SCRATCH_DIR"])/cytof/results/repfam/test-sim-6-8-3"
csvdir = "viz/csv"
mkpath(csvdir)

function get_best_lam(samples)
  lams = [s[:theta__lam] for s in samples]
  Ws = [s[:theta__W] for s in samples]
  Zs = [s[:theta__Z] for s in samples]
  I = length(lams[1])

  best_idx = [PlotUtils.estimate_ZWi_index(Zs, Ws, i) for i in 1:I]
  best_lam = [lams[best_idx[i]][i] for i in 1:I]

  return best_lam
end

for zind in 1:3
  for pmiss in (0.0, 0.6)
    for phi in (0.0, 100.0)
      # Get simulation name
      simname = "pmiss$(pmiss)-phi$(phi)-zind$(zind)"
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
      open("$(csvdir)/rfam-$(simname).csv", "w") do io
        writedlm(io, best_lam)
      end
    end
  end
end
