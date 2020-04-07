eye(n) = LinearAlgebra.I + zeros(n, n)

mutable struct PropCov
  value::Matrix{Float64}
  buffer::Matrix{Float64}
  iter::Int
  batchsize::Int
  scale::Float64
  momentum::Float64
  eps::Float64
end

function PropCov(value::Matrix{Float64}; iter::Int=0,
                 batchsize::Int=200, momentum::Float64=0.1,
                 eps::Float64=1e-6)
  # Get dimensions of prop cov
  K, K2 = size(value)
  @assert K == K2

  # Momentum needs to be between 0 and 1
  @assert 0 <= momentum <= 1

  # Create buffer
  buffer = zeros(Float64, batchsize, K)

  # Scale for proposal covariance
  scale = 2.38 * 2.38 / K

  return PropCov(value, buffer, iter, batchsize, scale, momentum, eps)
end

function PropCov(K::Int; sigma::Float64=1.0, iter::Int=0, batchsize::Int=200,
                 momentum::Float64=0.1, eps::Float64=1e-6)
  return PropCov(eye(K) * sigma, iter=iter, batchsize=batchsize,
                 momentum=momentum, eps=eps)
end

nparam(pc::PropCov) = size(pc.buffer, 2)

function update!(pc::PropCov, x::Vector{Float64})
  # Increment iteration
  pc.iter += 1

  # Compute current buffer index
  buffer_index = let
    idx = mod(pc.iter, pc.batchsize)
    idx == 0 ? pc.batchsize : idx
  end

  # Cache current value in buffer
  pc.buffer[buffer_index, :] .= deepcopy(x)

  if buffer_index == pc.batchsize
    # Update prop cov
    pc.value .= let
      S = cov(pc.buffer) * pc. scale * pc.momentum + pc.value * (1 - pc.momentum)
      S + eye(nparam(pc)) * pc.eps
    end
  end
end

# K = 30
# pc = Metropolis.PropCov(K, sigma=1e-2)
# for i in 1:10000 Metropolis.update!(pc, randn(K)) end
# pc.value
