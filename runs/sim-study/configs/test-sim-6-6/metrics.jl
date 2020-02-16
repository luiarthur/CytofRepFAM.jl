include(joinpath(@__DIR__, "../../../PlotUtils/PlotUtils.jl"))
include(joinpath(@__DIR__, "../../../PlotUtils/imports.jl"))

using Distributed

rmprocs(filter(w -> w > 1, workers()))
addprocs(27)

@everywhere include(joinpath(@__DIR__, "../../../PlotUtils/PlotUtils.jl"))
@everywhere include(joinpath(@__DIR__, "../../../PlotUtils/imports.jl"))

# Directory where results are
results_dir = "/scratchdata/alui2/cytof/results/repfam/test-sim-6-5-2"

### Create metrics ###

# Name of output file
@everywhere OUTPUT_FILE = "output.bson"

println("Producing metrics ...")

# directories to loop through
results_dirs = ["$(results_dir)/dataseed_$(dataseed)/mcmcseed_$(mcmcseed)/scale_$(scale)"
                for dataseed in 1:3
                for mcmcseed in 1:3
                for scale in (0, 1, 10)]

# Number of directories
num_dirs = length(results_dirs)

parsekey(T, key, path) = parse(T, match(Regex("(?<=$(key)_)\\d+"), path).match)

# Construct a data frame for R
function fill_metrics_df!(path, metrics, metrics_df)
  # Parse dataseed
  dataseed = parsekey(Int, "dataseed", path)

  # Parse mcmcseed
  mcmcseed = parsekey(Int, "mcmcseed", path)

  # Parse scale
  scale = parsekey(Int, "scale", path)

  I = length(metrics[:Rmean])
  for i in 1:I
    append!(metrics_df,
            DataFrames.DataFrame(i=i, 
                                 dataseed=dataseed,
                                 mcmcseed=mcmcseed,
                                 scale=scale,
                                 DIC=metrics[:DIC],
                                 LPML=metrics[:LPML],
                                 Rmean=metrics[:Rmean][i],
                                 R_975=metrics[:R_975][i],
                                 R_025=metrics[:R_025][i]))
  end
end

# Create make metrics function which only takes directory as input
@everywhere function getmetrics(d, OUTPUT_FILE)
  output = BSON.load("$(d)/$(OUTPUT_FILE)")
  metrics = Dict{Symbol, Any}()
  metrics[:LPML] = output[:LPML]
  metrics[:DIC] = output[:DIC]

  # Get posterior samples of W
  Ws = [o[:theta__W] for o in output[:samples][1]]

  # Posterior samples of R (I x K) -- number of active features / sample
  # Rs (I x num_mcmc_samples)
  Rs = hcat([vec(sum(W .> 0, dims=2)) for W in Ws]...)

  # Mean of R, by sample (vector of length I)
  Rmean = vec(mean(Rs, dims=2))

  # Mean of R, by sample (vector of length I)
  R_025 = PlotUtils.quantiles(Rs, .025, dims=2, drop=true)

  # Mean of R, by sample (vector of length I)
  R_975 = PlotUtils.quantiles(Rs, .975, dims=2, drop=true)

  # Append Rmean to metrics dictionary
  metrics[:Rmean] = Rmean

  # Append Rmean to metrics dictionary
  metrics[:R_025] = R_025

  # Append Rmean to metrics dictionary
  metrics[:R_975] = R_975

  return metrics
end

@everywhere make_metrics(d) = getmetrics(d, OUTPUT_FILE)

# Get all metrics in parallel
allmetrics = pmap(make_metrics, results_dirs)

# Put results in a data frame.
metrics_df = let
  df = DataFrames.DataFrame()
  for (d, metrics) in zip(results_dirs, allmetrics)
    fill_metrics_df!(d, metrics, df)
  end
  df
end

CSV.write("metrics.csv", metrics_df)

# Remove extra processors
rmprocs(filter(w -> w > 1, workers()))

println("DONE!")
