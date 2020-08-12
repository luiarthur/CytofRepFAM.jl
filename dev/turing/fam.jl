import Pkg; Pkg.activate("../repo/TuringBnpBenchmarks")

using Turing
using Distributions
using DelimitedFiles
import Random

idx, Y = let
  data = readdlm("patients.csv", ',')
  Int.(data[:, 1]), data[:, 2:end]
end

# FAM without missing data.
@model fam(y, idx, K, L0, L1, phi=0) = begin
  I = length(unique(idx))
  Nsum = length(idx)
  J = size(y, 2)

  alpha ~ Gamma(1, 1)
  v ~ filldist(Beta(alpha/K, 1), K)
  Z ~ arraydist([Bernoulli(v[k]) for j in 1:J, k in 1:K])
  lambda = tzeros(Int, Nsum)

  sigma ~ filldist(LogNormal(0, 1), I)
  delta0 ~ filldist(LogNormal(0, 0.1), L0)
  delta1 ~ filldist(LogNormal(0, 0.1), L1)
  mu0 = -cumsum(delta0)
  mu1 = cumsum(delta1)

  W ~ filldist(Dirichlet(K, 1/K), I)

  eta0 = Array{Float64, 3}(undef, I, J, L0)
  eta1 = Array{Float64, 3}(undef, I, J, L1)
  for i in 1:I
    for j in 1:J
      eta0[i, j, :] ~ Dirichlet(L0, 1)
      eta1[i, j, :] ~ Dirichlet(L1, 1)
    end
  end

  for j in 1:J
    for n in 1:Nsum
      i = idx[n] 

      lambda[n] ~ Categorical(W[:, i])
      z = Z[j, lambda[n]]
      Lz = (z == 0) ? L0 : L1
      muz = (z == 0) ? mu0 : mu1
      etaz = (z == 0) ? eta0 : eta1

      y[n, j] ~ UnivariateGMM(muz, fill(sigma[i], Lz), Categorical(etaz[i, j, :]))
    end
  end
end

Random.seed!(0)
@time chain = begin
    K = 25
    L0 = 5
    L1 = 5
    nburn = 10
    nsamps = 10
    iterations = nburn + nsamps
    sample(fam(Y, idx, K, L0, L1), Gibbs(PG(100)), iterations)
end;

nothing

#=
group(chain, :Z).value.data
=#
