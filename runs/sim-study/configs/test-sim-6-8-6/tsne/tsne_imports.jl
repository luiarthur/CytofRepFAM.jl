# ENV["PYCALL_JL_RUNTIME_PYTHON"] = "venv/bin/python"

import Pkg; Pkg.activate(joinpath(@__DIR__, "../../../../../"))
include(joinpath(@__DIR__, "../../../../PlotUtils/imports.jl"))
include(joinpath(@__DIR__, "../../../../PlotUtils/PlotUtils.jl"))
include(joinpath(@__DIR__, "../../../../PlotUtils/Population.jl"))
using DataFrames
using DelimitedFiles
using CSV
using PyCall
TSNE = PyCall.pyimport("sklearn.manifold").TSNE
const rfam = CytofRepFAM.Model

function impute(y; yQuantiles=[.0, .25, .5], pBounds=[.05, .8, .05], seed=1)
  I = length(y)

  # Missing mechanism parameters
  beta = hcat([rfam.gen_beta_est(vec(y[i]), yQuantiles, pBounds)
               for i in 1:I]...)

  # Mean of imputed values
  missMean = [beta[2, i] / (-2 * beta[1, i]) for i in 1:I]

  # Preimpute values
  y_imputed = deepcopy(y)
  for i in 1:I
    rfam.preimpute!(y_imputed[i], missMean[i])
  end

  return y_imputed
end

function get_best_lam(samples)
  n = length(samples) 
  Ws = [s[:theta__W] for s in samples]
  Zs = [s[:theta__Z] for s in samples]
  lams = [s[:theta__lam] for s in samples]
  I = length(lams[1])

  idx_best = [estimate_ZWi_index(Zs, Ws, i) for i in 1:I]
  return [lams[idx_best[i]][i] for i in 1:I]
end

function compute_tsne(path_to_simdat; use_complete_data=true, seed=0, verbose=2)
  simdat = BSON.load(path_to_simdat)[:simdat]

  if use_complete_data
    y = Matrix{Float64}.(simdat[:y_complete])
  else
    y = impute(simdat[:y])
  end

  # number of samples in y
  I = length(y)

  # Create TSNE objects
  tsne = [TSNE(verbose=verbose, random_state=seed) for _ in 1:I]

  # Fit TSNE for each sample
  foreach(i -> tsne[i].fit(y[i]), 1:I)

  return tsne
end

function compute_combined_tsne(path_to_simdat; use_complete_data=true, seed=0,
                               verbose=2, digits=5)
  # Load simulated data
  simdat = BSON.load(path_to_simdat)[:simdat]

  # Missing data indicator
  m = [isnan.(yi) for yi in simdat[:y]]

  # Determine what to do if missing data is present.
  if use_complete_data
    y = Matrix{Float64}.(simdat[:y_complete])
  else
    y = impute(simdat[:y])
  end

  # Missing indicator combined
  M = vcat(m...)

  # Combine samples into one matrix
  Y = vcat(y...)

  # number of samples in y
  I = length(y)

  # Sample sizes
  N = size.(y, 1)

  # Sample indicators
  sample_ind = vcat([fill(i, N[i]) for i in 1:I]...)

  # Create TSNE object
  tsne = TSNE(verbose=verbose, random_state=seed) 

  # Fit TSNE for samples combined
  tsne.fit(Y)

  # Get true labels
  true_labels = vcat(simdat[:lam]...)

  return (round.(tsne.embedding_, digits=digits), sample_ind,
          round.(Y, digits=digits), Int.(M), true_labels)
end

function gen_imgdir(path_to_output, suffix)
  d = splitdir(path_to_output)[1] 
  d1, d2 = splitdir(d)
  return joinpath(d1, "tsne", d2, suffix)
end
