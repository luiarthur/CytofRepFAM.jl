"""
Impute missing values in a column, using the other columns.

# Arguments
- Y: matrix (with preimputed values)
- M: missingness indictor (1 for missing, 0 for observed)
- k: column to impute
"""
function imputecolumn(Y, M, k)
  idx_mis = findall(M[:, k])
  idx_obs = findall(.!M[:, k])
  coef = Y[idx_obs, 1:end .!= k] \ Y[idx_obs, k]
  pred = Y[idx_mis, 1:end .!= k] * coef
  return pred
end


"""
Iteratively impute missing values in matrix, using the other columns.

# Arguments
- Y: matrix (with missing values)
- k: column to impute
"""
function impute(Y; maxiter=30, tol=1e-3, thresh=0, init=-3, verbose=1)
  @assert ndims(Y) == 2

  X = deepcopy(Y)
  M = isnan.(X)
  K = size(X, 2)
  
  # Initialize
  for k in 1:K
    num_missing = sum(M[:, k])
    neg_obs_idx = findall(X[:, k] .< 0)
    Xk_miss_idx = findall(M[:, k])
    X[Xk_miss_idx, k] .= init
  end

  # Difference of X between iterations
  diff_X = Float64[]

  for i in 1:maxiter
    X_old = deepcopy(X)
    if verbose
      print("\r$(i)/$(maxiter)")
    end

    for k in 1:K
      X[M[:, k], k] .= imputecolumn(X, M, k)
    end

    append!(diff_X, mean(abs.(X_old[M1] - X[M1])))

    if diff_X[end] < tol
      if verbose
        println("Convergence defected! Stopping early.")
      end
      break
    end
  end

  if verbose
    println()
  end

  return X, diff_X
end

# TEST
N, K = 1000, 3
Y = randn(N, K)
Y[Y .< -1] .= NaN
# @code_warntype imputecolumn(Y, M, 1)

# FIXME
impute(Y)
