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
  S = randn(P, P)
  SST = S * S'
  m = randn(P)
  (m, SST)
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

delta(n) = min(n^(-0.5), .01)
compute_acceptance_rate(x) = length(unique(x)) / length(x)

function adapt(init, epochs, batchsize; target_acceptance_rate=.23)
  x = copy(init)
  P = length(init)
  propCov = 0.01 * eye(P) / P
  for i in 1:epochs
    print("\r$(i)/$(epochs)")
    samps = gensamples(x, propCov, batchsize)
    acceptance_rate = compute_acceptance_rate(samps)
    # Adapt s, propCov
    propCov = let
      # See: http://probability.ca/jeff/ftpdir/adaptex.pdf
      Cest = cov(samps) # * (1 / i) + propCov * (i - 1) / i
      C1 = Cest * 5.66 / P + 1e-6 * eye(P)
      C0 = 0.01 * eye(P) / P
      if i <= 2 * P
        C = C0
      else
        C = 0.95 > rand() ? C1 : C0
      end
      C
    end
    x .= copy(samps[end])

    if i == epochs
      println()
      return (samps, propCov)
    end
  end
end


@time samps, propCov = adapt(randn(P), 1000, 500);
println("Dims: $P")
println("Accept rate: $(compute_acceptance_rate(samps))")
K = 5
println("Est. Cov.:")
pretty(cov(samps)[1:K, 1:K], digits=3)
println("True Cov.:")
pretty(cov(mvn)[1:K, 1:K], digits=3)
