# Generalized Bayes: https://arxiv.org/pdf/2006.05451.pdf

struct State
  cluster::Vector{Int}
  Z::Matrix{Bool} 
end

function loss(y::Matrix{Bool}, i::Integer, s::State)
  return sum(abs.(y[i, :] - s.Z[s.cluster[i]]))
end

function loss(y::Matrix{Bool}, s::State)
  N = length(s.cluster)
  return sum([loss(y, i, s) for i in 1:N])
end

function update_c!(s::State, y::Matrix{Bool})
  # TODO
end

function update_Z!(s::State, y::Matrix{Bool})
  # TODO
end
