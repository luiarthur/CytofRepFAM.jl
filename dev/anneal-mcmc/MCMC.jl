module MCMC

using Distributions
using LinearAlgebra

eye(n) = LinearAlgebra.I + zeros(n, n)

function _metropolis(x::Float64, logprob::Function, propsd::Float64)
  cand = curr + randn() * propsd
  met_ratio = logprob(cand) - logprob(curr)
  accept = met_ratio > log(rand())
  draw = accept ? cand : curr
  return (draw, accept)
end

function metropolis(x::Float64, logprob, propsd)
  return _metropolis(x, logprob, propsd)[1]
end

function metropolis(curr::Vector{Float64}, logprob::Function,
                    propcov::Matrix{Float64})
  cand = rand(MvNormal(curr, propcov))
  met_ratio = logprob(cand) - logprob(curr)
  out = met_ratio > log(rand()) ? cand : curr
  return out
end

function mcmc(init::Vector{T}, ll::Function, lp::Function;
              nmcmc::Int, nburn::Int, batchsize::Int=50,
              thin::Int=1,
              propcov_factor=1.0,
              propcov_init=eye(length(init)) * propcov_factor,
              discount=0.3,
              min_temper=1.0, max_temper=1.0) where T <: AbstractFloat
  @assert min_temper <= max_temper

  # Initialize MCMC
  state = deepcopy(init)

  # Number of parameters
  nparam = length(state)

  # Preallocate memory
  out = zeros(nmcmc, nparam)

  # Log density of target distribution, for a given temperature
  logprob(s::Vector{Float64}, temper::Float64) = ll(s) / temper + lp(s)

  # Total number of iterations
  total_iters = nmcmc + nburn

  # proposal covariance
  propcov = propcov_init

  # buffer
  buffer = zeros(batchsize, nparam)

  # Update function
  function update(s, iter, propcov)
    # Compute temperature
    power_temper = (nburn - iter + 1) / (nburn - 1)
    temper = iter < nburn ? max_temper^power_temper : min_temper

    # Log density of target distribution
    _logprob(s::Vector{Float64}) = logprob(s, temper)

    # Update state
    s .= metropolis(s, _logprob, propcov)
  end

  # Run MCMC
  for i in 1:total_iters
    print("\r $(i)/$(total_iters)")
    # Update state
    for t in 1:thin
      update(state, i, propcov)
    end

    # update buffer
    buffer_index = let
      idx = mod(i, batchsize)
      idx == 0 ? batchsize : idx
    end
    buffer[buffer_index, :] .= state

    # Update proposal covariance
    if buffer_index == batchsize
      # Update propcov
      propcov .= let
        a = cov(buffer) * discount 
        b = propcov * (1-discount) + eye(nparam) * 1e-6
        a + b
      end
    end

    # Keep samples
    if i > nburn
      out[i - nburn, :] .= state
    end
  end

  println()
  return out, propcov
end

end # module MCMC
