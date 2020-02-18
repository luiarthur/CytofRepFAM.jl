import Pkg; Pkg.activate(joinpath(@__DIR__, "../../"))  # CytofRepFam
include(joinpath(@__DIR__, "ParTemp.jl"))
using Distributions
using LinearAlgebra
import Random

function pretty(x; digits=nothing)
  if digits == nothing
    show(stdout, "text/plain", x)
  else
    show(stdout, "text/plain", round.(x, digits=digits))
  end
  println()
end

Random.seed!(0)

eye(n) = Matrix{Float64}(LinearAlgebra.I, n, n)

function gen_distribution(P)
  S = rand(InverseWishart(P, eye(P)))
  m = randn(P)
  (m, S)
end

# Dimensions
P = 10
m, S = gen_distribution(P)
mvn = MvNormal(S)


function gensamples(init::T, propCov, iters) where T
  logprob(x) = logpdf(mvn, x)
  out = T[]
  x = copy(init)
  for i in 1:iters
    x, _ = ParTemp.metropolisMV(x, logprob, propCov)
    append!(out, [copy(x)])
  end
  return out
end

delta(n) = min(n^(-0.5), .001)
compute_acceptance_rate(x) = length(unique(x)) / length(x)

function adapt(init, epochs, batchsize; target_acceptance_rate=.23, s=1.0)
  x = copy(init)
  propCov = s * eye(length(init))
  for i in 1:epochs
    print("\r$(i)/$(epochs)")
    samps = gensamples(x, propCov, batchsize)
    acceptance_rate = compute_acceptance_rate(samps)
    # Adapt s
    propCov = let
      factor = exp(delta(i))
      if acceptance_rate > target_acceptance_rate
        s *= factor
      else
        s /= factor
      end
      C = s * cov(hcat(samps...)') + 1e-6 * eye(P)
      C
    end
    x .= copy(samps[end])

    if i == epochs
      println()
      return (samps, propCov)
    end
  end
end


samps, propCov = adapt(randn(P), 500, 500, s=0.1);
println("Accept rate: $(compute_acceptance_rate(samps))")
K = 5
pretty(cov(samps)[1:K, 1:K])
pretty(cov(mvn)[1:K, 1:K], digits=3)
