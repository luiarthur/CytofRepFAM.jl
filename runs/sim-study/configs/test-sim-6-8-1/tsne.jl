# ENV["PYCALL_JL_RUNTIME_PYTHON"] = "venv/bin/python"
import Pkg; Pkg.activate(joinpath(@__DIR__, "../../../../"))

include(joinpath(@__DIR__, "../../../PlotUtils/imports.jl"))
using Pandas; const pd = Pandas

using PyCall
const TSNE = PyCall.pyimport("sklearn.manifold").TSNE

function get_best_lam(samples)
  n = length(samples) 
  Ws = [s[:theta__W] for s in samples]
  Zs = [s[:theta__Z] for s in samples]
  lams = [s[:theta__lam] for s in samples]
  I = length(lams[1])

  idx_best = [estimate_ZWi_index(Zs, Ws, i) for i in 1:I]
  return [Int64.(lams[idx_best[i]][i]) for i in 1:I]
end

# sns.pairplot(x_vars=:V1, y_vars=:V2, data=df, hue=:class,
#              plot_kws=edgecolor=Dict(:linewidth => 0, :s=> 13),
#              aspect=1, height=5);
# plt.savefig("test.pdf", bbox_inches="tight");
# plt.close();

# MAIN
results_dir = "/scratchdata/alui2/cytof/results/repfam/test-sim-6-8-1"
tsne_dir = joinpath(results_dir, "tsne")
path_to_simdat_1 = joinpath(results_dir, "pmiss0.6-phi0.0-zind1/simdat.bson")
path_to_simdat_2 = joinpath(results_dir, "pmiss0.6-phi0.0-zind2/simdat.bson")
path_to_simdat_3 = joinpath(results_dir, "pmiss0.6-phi0.0-zind3/simdat.bson")

phi = [0.0, 1.0, 10.0, 25.0]
zind = [1, 2, 3]

path_to_output_1 = joinpath(results_dir, "pmiss0.6-phi1.0-zind1/output.bson")
path_to_output_2 = joinpath(results_dir, "pmiss0.6-phi1.0-zind2/output.bson")
path_to_output_3 = joinpath(results_dir, "pmiss0.6-phi1.0-zind3/output.bson")

path_to_simdat = path_to_simdat_1 
path_to_output = path_to_output_1 
outdir = joinpath(tsne_dir, "z1")
use_complete_data = true
seed = 0
verbose = 2
markersize = 10
function compute_tsne(path_to_simdat, path_to_output, outdir;
                      use_complete_data=true, seed=0, verbose=2, 
                      markersize=10)
  simdat = BSON.load(path_to_simdat_1)[:simdat]
  output = BSON.load(path_to_output_1)

  mkpath(outdir)

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
    plt.savefig("$(outdir)/rfam$(i).pdf", bbox_inches="tight");
    plt.close();
  end
end


