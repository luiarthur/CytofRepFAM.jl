"""
Tuning parameter for metropolis step of a continuous parameter
with support on the real line.
"""
mutable struct TuningParam{T}
  value::T
  acceptanceCount::Int
  currentIter::Int
  batchsize::Int
  delta::Function
  targetAcc::Float64

  TuningParam(v::V) where V = new{V}(v, 0, 0, 50,
                                     n::Int -> max(n^(-0.5), 0.01),
                                     0.44)
end


"""
Performs update of the tuning parameter supplied, given whether
the proposed value was accepted or not.
"""
function update_tuning_param_default(param::TuningParam, accept::Bool)
  if accept
    param.acceptanceCount += 1
  end

  param.currentIter += 1

  if param.currentIter % param.batchsize == 0
    n = Int(floor(param.currentIter / param.batchsize))
    factor = exp(param.delta(n))
    if acceptanceRate(param) > param.targetAcc
      param.value *= factor
    else
      param.value /= factor
    end

    param.acceptanceCount = 0
  end

  return
end

function acceptanceRate(param::TuningParam)
  return param.acceptanceCount / param.batchsize
end

