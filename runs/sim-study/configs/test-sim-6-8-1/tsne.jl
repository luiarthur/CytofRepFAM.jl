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

function compute_tsne(path_to_simdat, path_to_output, imgdir;
                      method, use_complete_data=true, seed=0, verbose=2,
                      markersize=10)
  simdat = BSON.load(path_to_simdat)[:simdat]
  output = BSON.load(path_to_output)

  mkpath(imgdir)

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

  # Best lambdas
  best_lambda = get_best_lam(output[:samples][1])
  for i in 1:I
    tsne_df = pd.DataFrame(tsne[i].embedding_, columns=["comp1", "comp2"])
    tsne_df[Symbol("rfam")] = best_lambda[i]

    sns.pairplot(x_vars=:comp1, y_vars=:comp2, data=tsne_df, hue=:rfam,
                 plot_kws=edgecolor=Dict(:linewidth => 0, :s=> markersize),
                 aspect=1, height=5);
    plt.savefig("$(imgdir)/rfam$(i).pdf", bbox_inches="tight");
    plt.close();
  end
end


function gen_imgdir(path_to_output, suffix)
  d = splitdir(path_to_output)[1] 
  d1, d2 = splitdir(d)
  return joinpath(d1, "tsne", d2, suffix)
end

get_simdat_path(path_to_output) = joinpath(splitdir(path_to_output)[1], "simdat.bson")

## MAIN ###
results_dir = "/scratchdata/alui2/cytof/results/repfam/test-sim-6-8-1"

path_to_outputs = let
  path = [joinpath(results_dir, "pmiss$(pmiss)-phi$(phi)-zind$(zind)/output.bson")
          for phi in (0.0, 1.0),
              zind in (1, 2, 3),
              pmiss in (0.6)]  # pmiss in (0.0, 0.6)
  vec(path)
end

for path_to_output in path_to_outputs
  for method in ["rfam"]
    path_to_simdat = get_simdat_path(path_to_output)
    imgdir = gen_imgdir(path_to_output, method)
    mkpath(imgdir)

    println("Now on: ", imgdir)

    compute_tsne(path_to_simdat, path_to_output, imgdir;
                 use_complete_data=true, seed=0, verbose=2, 
                 markersize=10, method=method)
  end
end
