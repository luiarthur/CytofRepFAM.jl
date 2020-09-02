include("tsne_imports.jl")

get_output_path(path_to_output) = let
  joinpath(splitdir(path_to_output)[1], "output.bson")
end

# results_dir = "$(ENV["SCRATCH_DIR"])/cytof/results/repfam/test-sim-6-8-6"
results_dir = "../tex/results"
csvdir = "viz/csv"
mkpath(csvdir)

function get_best_lam(samples; true_Z=nothing, I=nothing)

  if typeof(samples) == String
    path_to_results = samples
    best_lam = [Int.(readdlm("$(path_to_results)/lam$(i)_best.txt")) 
                for i in 1:I]
    best_Z = [Int.(readdlm("$(path_to_results)/Z$(i)_best.txt"))
              for i in 1:I]
  else
    lams = [s[:theta__lam] for s in samples]
    Ws = [s[:theta__W] for s in samples]
    Zs = [s[:theta__Z] for s in samples]
    I = length(lams[1])

    best_idx = [PlotUtils.estimate_ZWi_index(Zs, Ws, i) for i in 1:I]
    best_lam = [lams[best_idx[i]][i] for i in 1:I]
    best_Z = [Zs[best_idx[i]] for i in 1:I]
  end

  if true_Z == nothing
    return best_lam
  else
    population = Population()
    foreach(z -> label(population, z), eachcol(true_Z))
    relabeled_lam = [[label(population, best_Z[i][:, lam_in])
                      for lam_in in best_lam[i]] for i in 1:I]

    return relabeled_lam
  end
end

for zind in 1:3
  for pmiss in (0.0, 0.2)
    for phi in (0, 1, 10, 25, 100)
      # Get simulation name
      simname = "pmiss$(pmiss)-phi$(phi)-zind$(zind)"
      println("Processing: ", simname)

      # Get path to simuldate data.
      # path_to_output = joinpath(results_dir, simname, "output.bson")

      # Get output 
      # output = BSON.load(path_to_output)
      
      # Get relevant parameters
      # samples =  [s[1] for s in output[:samples][1]]

      # Get ture Z
      true_Z = readdlm("../Z$(zind).txt", comments=true, comment_char='#')

      # Get best lam (for each sample)
      # best_lam = vcat(get_best_lam(samples, true_Z=true_Z)...)
      path_to_results = joinpath(results_dir, simname, "img/txt")
      best_lam = vcat(get_best_lam(path_to_results, true_Z=true_Z, I=2)...)

      # Write output to file.
      open("$(csvdir)/rfam-$(simname).csv", "w") do io
        writedlm(io, best_lam)
      end
    end
  end
end
