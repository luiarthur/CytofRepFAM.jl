# ENV["PYCALL_JL_RUNTIME_PYTHON"] = "venv/bin/python"
import Pkg; Pkg.activate(joinpath(@__DIR__, "../../../../"))

include(joinpath(@__DIR__, "../../../PlotUtils/imports.jl"))
using Pandas; const pd = Pandas

using PyCall
TSNE = PyCall.pyimport("sklearn.manifold").TSNE

function get_best_lam(samples)
  n = length(samples) 
  Ws = [s[:theta__W] for s in samples]
  Zs = [s[:theta__Z] for s in samples]
  lams = [s[:theta__lam] for s in samples]
  I = length(lams[1])

  idx_best = [estimate_ZWi_index(Zs, Ws, i) for i in 1:I]
  return [Int64.(lams[idx_best[i]][i]) for i in 1:I]
end

function compute_tsne(path_to_simdat; use_complete_data=true, seed=0, verbose=2)
  simdat = BSON.load(path_to_simdat)[:simdat]

  if use_complete_data
    y = Matrix{Float64}.(simdat[:y_complete])
  else
    # TODO
  end

  # number of samples in y
  I = length(y)

  # Create TSNE objects
  tsne = [TSNE(verbose=verbose, random_state=seed) for _ in 1:I]

  # Fit TSNE for each sample
  foreach(i -> tsne[i].fit(y[i]), 1:I)

  return tsne
end

function plot_tsne(tsne, path_to_output, imgdir; method, markersize=10)
  output = BSON.load(path_to_output)
  mkpath(imgdir)

  # number of samples in y
  I = length(tsne)

  # Best lambdas
  best_lambda = get_best_lam(output[:samples][1])

  # Noisy class couunts (shouold be 0)
  num_noisy = zeros(Int, I)

  for i in 1:I
    # This shouold be zero
    num_noisy[i] = sum(best_lambda[i] .== 0)

    tsne_df = pd.DataFrame(tsne[i].embedding_, columns=["comp1", "comp2"])
    tsne_df[Symbol("rfam")] = best_lambda[i]

    sns.pairplot(x_vars=:comp1, y_vars=:comp2, data=tsne_df, hue=:rfam,
                 plot_kws=edgecolor=Dict(:linewidth => 0, :s=> markersize),
                 aspect=1, height=5);
    plt.savefig("$(imgdir)/rfam$(i).pdf", bbox_inches="tight");
    plt.close();
  end

  open("$(imgdir)/num_noisy.txt", "w") do io
    println(io, num_noisy)
  end
end


function gen_imgdir(path_to_output, suffix)
  d = splitdir(path_to_output)[1] 
  d1, d2 = splitdir(d)
  return joinpath(d1, "tsne", d2, suffix)
end


## MAIN ###
get_simdat_path(path_to_output) = joinpath(splitdir(path_to_output)[1], "simdat.bson")
results_dir = "/scratchdata/alui2/cytof/results/repfam/test-sim-6-8-1"

tsne = [let
          path = joinpath(results_dir, "pmiss0.0-phi0.0-zind$(zind)/simdat.bson")
          compute_tsne(path, use_complete_data=true, seed=0, verbose=2)
        end for zind in (1, 2, 3)]

path_to_outputs = let
  path = [joinpath(results_dir, "pmiss$(pmiss)-phi$(phi)-zind$(zind)/output.bson")
          for phi in (0.0, 1.0, 10.0, 25.0),
              zind in (1, 2, 3),
              pmiss in (0.6, 0.0)]  # pmiss in (0.0, 0.6)
  vec(path)
end

for path_to_output in path_to_outputs
  for method in ["rfam"]
    imgdir = gen_imgdir(path_to_output, method)
    mkpath(imgdir)

    println("Now on: ", imgdir)

    zind = let
      m = match(r"(?<=zind)\d+", path_to_output).match
      parse(Int, m)
    end

    plot_tsne(tsne[zind], path_to_output, imgdir;
              markersize=10, method=method)
  end
end
