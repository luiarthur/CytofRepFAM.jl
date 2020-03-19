"""
m: minibatch size
N: full data size
"""
function subsample_index(m::Integer, N::Integer)
  @assert m < N
  shuffled_idx = Random.shuffle(1:N)
  return shuffled_idx[1:m], shuffled_idx[(m+1):end]
end

function sample_minibatch(batchsizes::Vector{S},
                          y::Vector{Matrix{T}}) where {T <: AbstractFloat,
                                                       S <: Integer}
  I = length(y)
  @assert length(batchsizes) == I
  N = size.(y, 1)
  J = size(y[1], 2)

  idx_mini, idx_mega = zip([subsample_index(batchsizes[i], N[i]) for i in 1:I]...)
  ymini = [y[i][idx_mini[i], :] for i in 1:I]
  ymega = [y[i][idx_mega[i], :] for i in 1:I]

  return idx_mini, idx_mega, Data(ymini), Data(ymega)
end
