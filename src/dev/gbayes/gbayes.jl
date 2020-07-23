# Generalized Bayes: https://arxiv.org/pdf/2006.05451.pdf

using Distributions

struct State
  cluster::Vector{Int}
  Z::Matrix{Bool} 
  v::Vector{Float64}
  alpha::Float64
end

function loss(yi::Vector{Bool}, zk::Vector{bool})
  return sum(abs.(yi - zk))
end

function loss(y::Matrix{Bool}, i::Integer, s::State)
  return loss(y[i, :], s.Z[s.cluster[i]])
end

function loss(y::Matrix{Bool}, s::State)
  N = length(s.cluster)
  return sum([loss(y, i, s) for i in 1:N])
end

function log_prior_Z(s::State, phi::Float64)
  # TODO
end

function update_c!(s::State, i::Integer, y::Matrix{Bool}, lam::Float64)
  K = size(s.Z, 2)
  log_prob = [-lam * loss(y[i], s.Z[:, k]) for k in 1:K]
  unnormalized_prob = exp.(log_prob - max(log_prob))
  prob = unnormalized_prob ./ sum(unnormalized_prob)
  s.cluster[i] = wsample(prob)
  return
end

function update_c!(s::State, y::Matrix{Bool}, lam::Float64)
  N = length(s.cluster)
  for i in 1:N
    update_c!(s, i, y, lam)
  end
end

function update_Z!(s::State, j:Integer, k::Integer, y::Matrix{Bool}, lam::Float64)
  # TODO
end

function update_Z!(s::State, y::Matrix{Bool}, lam::Float64)
  J, K = size(s.Z)
  for j in 1:J
    for k in 1:K
      update_Z(s, j, k, y, lam)
    end
  end
end
