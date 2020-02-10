module Regress
using StatsBase

"""
Impute missing values in a column, using the other columns.

# Arguments
- Y: matrix (with preimputed values)
- M: missingness indictor (1 for missing, 0 for observed)
- k: column to impute
- thresh: In regression, only train model with responses less than this value.
"""
function imputecolumn(Y, M, k; thresh=Inf)
  idx_mis = findall(M[:, k])
  idx_obs = findall((.!M[:, k]) .& (Y[:, k] .< thresh))
  coef = Y[idx_obs, 1:end .!= k] \ Y[idx_obs, k]
  pred = Y[idx_mis, 1:end .!= k] * coef
  return pred
end


"""
Iteratively impute missing values in matrix, using the other columns.

# Arguments
- Y: matrix (with missing values)
- maxiter: Maximum number of iterations to impute missing values.
- tol: if difference between current and next matrices are less than this,
  declare convergence!
- init: initialize missing values at this.
- verbose: if value > 0, print some information about the algorithm.
- thresh: See `imputecolumn`.
"""
function impute(Y; maxiter=30, tol=1e-3, init=-3, verbose=1, thresh=Inf)
  @assert ndims(Y) == 2

  X = deepcopy(Y)
  M = isnan.(X)
  K = size(X, 2)
  
  # Initialize
  X[M] .= init

  # Difference of X between iterations
  diff_X = Float64[]

  for i in 1:maxiter
    X_old = deepcopy(X)
    if verbose > 0
      print("\r$(i)/$(maxiter)")
    end

    for k in 1:K
      X[M[:, k], k] .= imputecolumn(X, M, k, thresh=thresh)
    end

    append!(diff_X, mean(abs.(X_old[M] - X[M])))

    if diff_X[end] < tol
      if verbose > 0
        println()
        println("Convergence defected! Stopping early.")
      end
      break
    end
  end

  if verbose > 0
    println()
  end

  return X, diff_X
end

end  # module Regress

#= Demo
using CSV
using Statistics
using DataFrames
using PyPlot; const plt = PyPlot.plt
Y = coalesce.(CSV.read("../../runs/data/cb_transformed.csv"), NaN)
Y_mat = Matrix(Y[Y[!, :sample_id] .== 1, Not(:sample_id)])

Y_complete, diff_Y = Regress.impute(Y_mat, init=-6, thresh=.5, maxiter=50)

function plotme(j; magthresh=7)
  yj = Y_complete[isnan.(Y_mat[:, j]), j]
  plt.hist(yj[abs.(yj) .< magthresh , :], bins=50, alpha=.5)
  axvline(0, lw=2)
  println("Proportion of missing: $(mean(isnan.(Y_mat[:, j])))")
end

plotme(2)  # 16 is interesting
plt.close();
=#

