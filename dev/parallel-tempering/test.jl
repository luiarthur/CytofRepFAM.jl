# using RCall
# @rimport rcommon

include("imports.jl")
using Distributed
rmprocs(filter(w -> w > 1, workers()))
addprocs(16)

@everywhere include("imports.jl")
@everywhere Random.seed!(0)

@everywhere eye(T::Type, n::Int) = Matrix{T}(LinearAlgebra.I, n, n)

@everywhere function gen_distribution(P)
  S = rand(InverseWishart(P, eye(Float64, P)))
  m = randn(P) * 10
  (m, S)
end

# True distribution has this covariance
# @everywhere P = 5  # NOTE: works great
@everywhere P = 20 # NOTE: works
# @everywhere P = 50 # NOTE: dimension too high?
@everywhere m, S = gen_distribution(P)
@everywhere mvn = MvNormal(m, S)

@everywhere function update(init::T, tau, propCov, iters) where T
  samps = T[]
  append!(samps, [copy(init)])

  logprob(z) = logpdf(mvn, z) / tau

  x = copy(init)
  for _ in 1:iters
    x = ParTemp.metropolisMV(x, logprob, propCov)[1]
    append!(samps, [copy(x)])
  end

  return samps
end

taus = 1.1 .^ (0:15)
# taus = .99 .^ (15:-1:0)
iters = 100
epochs = 100
buffsize = 1000


@time samps = let
  out = nothing
  state = [rand(P) for _ in taus]
  propCov = [1e-6 * eye(Float64, P) for _ in taus]
  states = [ParTemp.Buffer(Vector{Float64}, buffsize) for _ in taus]

  for i in 1:epochs
    print("\rIter $i / $epochs")
    out = pmap(update,
               state,
               taus,
               propCov,
               [iters for _ in taus])

    for t in 1:length(taus)
      ParTemp.append!(states[t], out[t])
    end

    propCov = [cov(states[t].x) + 1e-6 * eye(Float64, P)
               for t in 1:length(taus)]
    state = map(o -> o[end], out)

    # Swap chains
    for t in length(taus):-1:2
      j = t
      i = j - 1
      logprob(k) = logpdf(mvn, state[k])
      logprobswap = (1/taus[j] - 1/taus[i]) * (logprob(i) - logprob(j))
      if logprobswap > log(rand())
        println("Swapping ($j, $i)")
        tmp = copy(state[j])
        state[j] = copy(state[i])
        state[i] = tmp
      end
    end
  end
  println()

  states[1].x
end

begin
  tail = Int(length(samps) / 2)
  K = 5
  mean(samps[(end-tail):end])
  println("Est. Cov:")
  show(stdout, "text/plain", cov(samps[(end-tail):end])[1:K, 1:K])
  # show(stdout, "text/plain", cor(hcat(samps[(end-tail):end]...)')[1:K, 1:K])

  println()
  println()

  println("True Cov:")
  show(stdout, "text/plain", cov(mvn)[1:K, 1:K])
  # show(stdout, "text/plain", cor(mvn)[1:K, 1:K])
  # rcommon.plotPosts(Matrix(hcat(samps[(end-tail):end]...)')[:, 1:K]);

  println("mean(mvn), mean(samps)")
  show(stdout, "text/plain", [mean(mvn) mean(samps)])

  nothing
end;
