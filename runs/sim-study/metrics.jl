include(joinpath(@__DIR__, "../PlotUtils/PlotUtils.jl"))
include(joinpath(@__DIR__, "../PlotUtils/imports.jl"))

using Distributed

rmprocs(filter(w -> w > 1, workers()))
addprocs(15)

@everywhere include(joinpath(@__DIR__, "../PlotUtils/PlotUtils.jl"))
@everywhere include(joinpath(@__DIR__, "../PlotUtils/imports.jl"))

if length(ARGS) == 0
  results_dir = "/scratchdata/alui2/cytof/results/repfam/test-sim-6-5"
else
  results_dir = ARGS[1]  # path to results directory
end

### Create metrics ###

# Name of output file
@everywhere OUTPUT_FILE = "output.bson"

println("Producing metrics ...")

# directories to loop through
results_dirs = ["$(results_dir)/seed_$(seed)/scale_$(scale)"
                for seed in 1:3
                for scale in (0, 1, 10)]

# Number of directories
num_dirs = length(results_dirs)

parsekey(T, key, path) = parse(T, match(Regex("(?<=$(key)_)\\d+"), path).match)

# Construct a data frame for R
function fill_R_df!(path, metrics, R_df)
  # Parse seed
  seed = parsekey(Int, "seed", path)

  # Parse scale
  scale = parsekey(Int, "scale", path)

  for metric in metrics
    Kmcmc, m = metric
    I = length(m[:Rmean])
    for i in 1:I
      append!(R_df, DataFrames.DataFrame(Kmcmc=Kmcmc, i=i,
                                         seed=seed, scale=scale,
                                         Rmean=m[:Rmean][i],
                                         R_975=m[:R_975][i],
                                         R_025=m[:R_025][i]))
    end
  end
end

# Create make metrics function which only takes directory as input
@everywhere make_metrics(d) = PlotUtils.make_metrics(d, OUTPUT_FILE, thresh=.01)

# Get all metrics in parallel
allmetrics = pmap(make_metrics, results_dirs)

# Put results in a data frame.
R_df = let
  _R_df = DataFrames.DataFrame()
  for (d, metrics) in zip(results_dirs, allmetrics)
    fill_R_df!(d, metrics, _R_df)
  end
  _R_df
end

# TODO: Plot R-metrics nicely
# for each seed, make the R by Kmcmc graphs
for R_df_seed in groupby(R_df, :seed)
  seed = unique(R_df_seed.seed)[1]
  unique_i = unique(R_df_seed.i)
  I = length(unique_i)
  plt.figure()
  for i in 1:I
    R_df_seed_i = R_df_seed[R_df_seed.i .== i, :]
    plt.subplot(I, 1, i)
    # Put graphs for different scales in same graph for same seed
    for R_df_seed_i_scale in groupby(R_df_seed_i, :scale)
      scale = unique(R_df_seed_i_scale.scale)[1]
      df = sort(R_df_seed_i_scale, :Kmcmc)
      plt.plot(df.Kmcmc, df.Rmean, label="scale: $(scale)", lw=2, marker="o")
      plt.fill_between(df.Kmcmc, df.R_025, df.R_975, alpha=.5)
    end
    plt.ylabel("R$(i)")
  end

  outdir = "$(results_dir)/metrics/seed_$(seed)"
  mkpath(outdir)

  # plt.legend(loc="lower right")
  plt.legend(loc="best")
  plt.xlabel("K")
  plt.savefig("$(outdir)/R.pdf", bbox_inches="tight")
  plt.close()
end

# Remove extra processors
rmprocs(filter(w -> w > 1, workers()))

println("DONE!")
