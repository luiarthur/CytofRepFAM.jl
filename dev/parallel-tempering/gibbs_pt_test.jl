# NOTE: Currently, not used.
module TestGibbsPT

import Pkg; Pkg.activate("../..")  # CytofRepFAM
import CytofRepFAM.MCMC, Random
using Distributions, Distributed

mutable struct Param
  mu
  sig2
  w
  lam
end


struct Priors
  mu
  sig2
  w
end

Priors(K) = Priors(Normal(0, 1), InverseGamma(3, 2), Dirichlet(K, 1/K))

Param(N, priors) = let
  K = length(priors.w.alpha)
  Param(rand(priors.mu, K), rand(priors.sig2),
        rand(priors.w), rand(1:K, N))
end

# Set up number of cores to use
maxcores = 20

# if length(workers()) != maxcores
#   rmprocs(filter(w -> w > 1, workers()))
#   addprocs(maxcores)
#   @everywhere begin
#     import Pkg; Pkg.activate("../")
#     using CytofRepFAM
#   end
# end

function update_mu(k, theta, priors, y)
  idx_k = findall(theta.lam .== k)
  nk = length(idx_k)
  yk = y[idx_k]
  newPrec = 1/var(priors.mu) + nk/theta.sig2
  newVar = 1/newPrec
  newMean = newVar * (mean(priors.mu)/var(priors.mu) + sum(yk)/theta.sig2)
  theta.mu[k] = rand(Normal(newMean, sqrt(newVar)))
end

function update_mu(theta, priors, y)
  for k in 1:length(theta.mu)
    update_mu(k, theta, priors, y)
  end
end

function update_sig2(theta, priors, y)
  N = length(y)
  ss = sum((y - theta.mu[theta.lam]) .^ 2)
  newShape = N / 2 + shape(priors.sig2)
  newScale = ss / 2 + scale(priors.sig2)
  theta.sig2 = rand(InverseGamma(newShape, newScale))
end

function update_w(theta, priors, y)
  N = length(y)
  K = length(theta.w)
  counts = zeros(Int, K)
  for i in 1:N
    counts[theta.lam[i]] += 1
  end
  newAlpha = priors.w.alpha + counts
  theta.w = rand(Dirichlet(newAlpha))
end

function update_lam(theta, priors, y)
  N = length(y)
  for i in 1:N
    update_lam(i, theta, priors, y)
  end
end

function update_lam(i, theta, priors, y)
  lls = MCMC.lpdf_normal.(y[i], theta.mu, sqrt(theta.sig2))
  lp = log.(theta.w)
  logprob = lls + lp
  theta.lam[i] = MCMC.wsample_logprob(logprob)
end

function update(theta, priors, y)
  update_mu(theta, priors, y)
  update_sig2(theta, priors, y)
  update_lam(theta, priors, y)
  update_w(theta, priors, y)
end


begin
  mu_true = [-2, 3, 1]
  sig2_true = 0.01
  w_true = [.3, .5, .2]
  K_true = length(mu_true)
  nobs = 100
  mus = wsample(mu_true, w_true, nobs)
  y = rand.(Normal.(mus, sqrt(sig2_true)))

  K = 5
  priors = Priors(5)
  theta = Param(nobs, priors)

  nburn = 100
  nmcmc = 1000
  states = Param[]

  for i in 1:nburn
    update(theta, priors, y)
  end

  for i in 1:nmcmc
    update(theta, priors, y)
    append!(states, [deepcopy(theta)])
  end

  extract(sym) = [getfield(s, sym) for s in states]

  mus = extract(:mu)
  sig2s = extract(:sig2)
  ws = extract(:w)

  println("mu: ", mean(mus))
  println("sig2: ", mean(sig2s))
  println("w: ", mean(ws))
end

end
