module ImportanceSampling

using Distributions
using LinearAlgebra: norm

function pairwise_dist(Z)
  K = size(Z, 2)
  return [norm(Z[:, k] - Z[:, h], 1) for k in 2:K for h in 1:(k - 1)]
end

function rep_fn(Z, phi)
  d = pairwise_dist(Z)
  if phi == 0
    return 1.0
  else
    log_rep_term = phi * sum(log1p.(-exp.(-d)))
    rep_term = exp(log_rep_term)
    return rep_term
  end
end

"""
Approximate the expectation of a scalar function of Z ~ rFAM E[f(Z)|phi] using
importance sampling.
"""
function exp_f_rfam(; N::Integer, J::Integer, K::Integer,
                    repulsive_fn::Function, f::Function)
  # NOTE: In the importance sampling, the density of the regular
  # IBP component is cancelled, leaving the repulsive component.
  # Therefore, we only need binary random matrices to approximate
  # mean functions of Z ~ rfam.
  Zs = [rand(Bool, J, K) for _ in 1:N]

  weights = repulsive_fn.(Zs)
  numer = sum(f.(Zs) .* weights)
  denom = sum(weights)

  if denom == 0
    println("WARNING: weights are all 0!")
    return 0.0
  else
    return numer / denom
  end
end

min_pairwise_dist_ge_u(Z, u) = minimum(pairwise_dist(Z)) .>= u

end  # ImportanceSampling
