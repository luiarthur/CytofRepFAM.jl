#=
import Pkg; Pkg.activate(".")
=# 
using Distributions
import Random
using RCall
import LinearAlgebra

include("MCMC.jl")

prettymat(X) = begin show(stdout, "text/plain", X); println(); end

### MAIN ###
Random.seed!(1)

# easy = true  # bivariate normal
easy = false  # 250-variate normal

if easy 
  K = 3
  S = MCMC.eye(K)
  S[1,2] = S[2,1] = -0.8
  S[1,3] = S[3,1] = 0.5
else
  K = 50
  S = rand(InverseWishart(K, MCMC.eye(K)))
end

mvn = MvNormal(randn(K) * 10, S)
logprob(x::Vector{Float64}) = logpdf(mvn, x)

init = randn(K) 

lp(s::Vector{Float64}) = 0.0
out, propcov = MCMC.mcmc(init, logprob, lp, thin=10, batchsize=100, 
                         temper_power=3,
                         propcov_init=MCMC.eye(length(init)) * 1e-3,
                         nburn=10000, nmcmc=1000, max_temper=1000.0,
                         verbose=1)

acc_rate = size(unique(out, dims=1), 1) / size(out, 1)

ll = [logprob(out[i, :]) for i in 1:size(out, 1)]
R"plot($ll, type='l')"

G = 8
println("cov(out):")
prettymat(round.(cov(out[:, 1:G]), digits=3))
println("cov(mvn):")
prettymat(round.(cov(mvn)[1:G, 1:G], digits=3))

[mean(out, dims=1)' mean(mvn)]

R"rcommon::plotPosts($(out)[, 1:5])";
println("Acceptance rate: $(acc_rate)")
