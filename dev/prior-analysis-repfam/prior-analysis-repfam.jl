import Pkg; Pkg.activate("../../")
include("HMC.jl")
using Distributions
using StatsBase

function grad_logprobH(H, j, k, phi)
  J, K = size(H)
  out = -H[j, k]
  for m in 1:K
    if m != k
      dmk = mean((H[:, k] - H[:, m]) .^ 2)
      tmp = zero(H[1])
      tmp += mean(H[:, k] - H[:, m]) / (exp(dmk / phi) - 1)
      out += tmp * 2 / phi
    end
  end
end

function grad_logprobH(H, phi)
  J, K = size(H)
  return [grad_logprobH(H, j, k, phi) for j in 1:J, k in 1:K]
end

"""
log1m(x) = log1p(-x) = log(1 - x) 
"""
log1m(x) = log1p(-x)

function logprobH(H, phi)
  J, K = size(H)
  normal_kernel = -sum((H .^ 2)) / 2
  penalty = zero(H[1])

  for m in 1:(K - 1)
    for n in (m + 1):K
      d = mean((H[:, m] - H[:, n]) .^ 2)
      penalty += log1m(exp(-d / phi))
    end
  end

  return normal_kernel + penalty
end

struct State H end

function drawH(J, K, phi)
  if phi == 0
    randn(J, K)
  else
    # TODO
    0
  end
end

function drawH(n, J, K, phi)
  if phi == 0
    [randn(J, K) for _ in 1:n]
  else
    # TODO
    
  end
end


"""
NOTE:
Let H be a JxK continuous matrix, with H_{jk} \\in \\mathbb{R}.
Let Z be a JxK binary matrix.
Let v be a K-dimensional continuous vector, bounded within the unit interval.

Model:
Z_{jk} = 1(v_{k} > \\Phi(H_{jk}))
v_k | \\alpha ~ Beta(\\alpha/K, 1)
p(H) \\propto \\prod_{j,k} exp(-H_{jk}^2/2) * 
              \\prod{k\\ne k'} 1 - \\exp(-d_{k,k'} / \\phi)
"""
function drawZ(J, K, alpha, phi; stickbreak=false)
  if stickbreak
    v = cumprod(rand(Beta(alpha, 1), 1, K), dims=2)
  else
    v = rand(Beta(alpha / K, 1), 1, K)
  end
  H = drawH(J, K, phi)
  return v .> cdf.(Normal(0, 1), H)
end

# Simulate
ndraws = 100000
Zs = [drawZ(2, 2, 1.0, 0.0, stickbreak=true) for _ in 1:ndraws]
countmap(Zs)
