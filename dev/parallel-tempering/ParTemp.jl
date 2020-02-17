module ParTemp

import Base.append!
using Distributions

function metropolisMV(curr, logprob, S)
  cand = rand(MvNormal(curr, S))
  logU = log(rand())
  logP = logprob(cand) - logprob(curr)
  accept = logP > logU
  draw = accept ? cand : curr
  return (draw, accept)
end

mutable struct Buffer{T}
  n::Int
  x::Vector{T}
  Buffer(T::Type, n::Int) = new{T}(n, T[])
end

tail(x, n) = n > length(x) ? x : x[(end-n+1):end]

function append!(buffer::Buffer, x)
  append!(buffer.x, x)
  buffer.x = tail(buffer.x, buffer.n)
end

end
