# NOTE
# Did a prior analysis for a version of repFAM where Z is
# a function of H, v. The results were underwhelming, and we couldn't
# achieve repulsiveness. Need to give up this idea.

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
  return out
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

function drawH(n, J, K, phi; num_leapfrog_steps=2^6, eps=.01, burn=0)
  lp(H) = logprobH(H, phi)
  grad_lp(H) = grad_logprobH(H, phi)

  H = randn(J, K)
  for _ in 1:burn
    H, _ = HMC.hmc_update(H, lp, grad_lp, num_leapfrog_steps, eps)
  end

  Hs = [zeros(J, K) for _ in 1:n]
  lpHs = zeros(n)
  
  for i in 1:n
    H, lpH = HMC.hmc_update(H, lp, grad_lp, num_leapfrog_steps, eps)
    Hs[i] .= copy(H)
    lpHs[i] = lpH
  end

  return Hs, lpHs
end

function drawZ(N, J, K, alpha, phi;
               stickbreak=false, num_leapfrog_steps=2^6, eps=.01, burn=0)
  if stickbreak
    vs = [cumprod(rand(Beta(alpha, 1), 1, K), dims=2) for _ in 1:N]
  else
    vs = [rand(Beta(alpha / K, 1), 1, K) for _ in 1:N]
  end

  Hs, lpH = drawH(N, J, K, phi,
                  num_leapfrog_steps=num_leapfrog_steps, eps=eps, burn=burn)

  return [vs[i] .> cdf.(Normal(0, 1), Hs[i]) for i in 1:N]
end


# Simulate
N, alpha, phi, J, K, burn = (10000, 0.5, 100.0, 2, 2, 100)
@time Zs = drawZ(N, J, K, alpha, phi, stickbreak=true, burn=burn)
Zcounts = countmap(Zs)
sort(collect(Zcounts), by=x -> -x.second)
