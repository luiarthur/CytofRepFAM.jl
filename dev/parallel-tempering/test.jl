include("ParTemp.jl")
# using RCall
# @rimport rcommon

using Distributions
using LinearAlgebra
import Random

Random.seed!(2)

function gen_distribution(P)
  rand(InverseWishart(P, eye(Float64, P)))
end

# True distribution has this covariance
# P = 5  # NOTE: works great!
P = 50   # FIXME: Doesn't work for high dims?
S = gen_distribution(P)

eye(T::Type, n::Int) = Matrix{T}(LinearAlgebra.I, n, n)


function update(init::T, tau, propCov, iters) where T
  samps = T[]
  append!(samps, [copy(init)])

  logprob(z) = logpdf(MvNormal(S), z) / tau

  x = copy(init)
  for _ in 1:iters
    x = ParTemp.metropolisMV(x, logprob, propCov)[1]
    append!(samps, [copy(x)])
  end

  return samps
end

taus = 1.5 .^ (0:15)
iters = 50
epochs = 200

@time samps = let
  out = nothing
  state = [rand(P) for _ in taus]
  propCov = [1e-6 * eye(Float64, P) for _ in taus]
  state1 = Vector{Float64}[]
  for i in 1:epochs
    print("\rIter $i / $epochs")
    out = map(update,
              state,
              taus,
              propCov,
              [iters for _ in taus])
    # FIXME. What's a good way to do this?
    # How to estimate a good proposal?
    propCov = [cov(out[t]) * 1/i + propCov[t] * (i - 1)/i +
               1e-6 * eye(Float64, P)
               for t in 1:length(taus)]
    state = map(o -> o[end], out)

    # Swap chains
    for t in length(taus):-1:2
      j = t
      i = j -1
      logprob(k) = logpdf(MvNormal(S), state[k])
      probswap = (1/taus[j] - 1/taus[i]) * (logprob(i) - logprob(j))
      if probswap > log(rand())
        println("Swapping ($j, $i)")
        tmp = copy(state[j])
        state[j] = copy(state[i])
        state[i] = tmp

        # tmp = propCov[j]
        # propCov[j] = propCov[i]
        # propCov[i] = tmp
      end
    end

    append!(state1, out[1])
    # propCov[1] = let
    #   # S = hcat(state1...)'
    #   # cor(S + randn(size(S)) * 1e-6)
    #   cov(state1) + eye(Float64, P) * 1e-6
    # end
  end

  println()
  state1
end

begin
  tail = Int(length(samps) / 2)
  K = 3
  mean(samps[(end-tail):end])
  println("Est. Cov:")
  show(stdout, "text/plain", cov(samps[(end-tail):end])[1:K, 1:K])

  println()
  println()

  println("True Cov:")
  show(stdout, "text/plain", S[1:K, 1:K])
  rcommon.plotPosts(Matrix(hcat(samps[(end-tail):end]...)')[:, 1:K]);

  nothing
end;
