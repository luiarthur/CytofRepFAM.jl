module MCMC

using Distributions
using LinearAlgebra
include("PropCov.jl")

function _metropolis(x::Float64, logprob::Function, propsd::Float64)
  cand = curr + randn() * propsd
  met_ratio = logprob(cand) - logprob(curr)
  accept = met_ratio > log(rand())
  draw = accept ? cand : curr
  return (draw, accept)
end

function metropolis(x::Float64, logprob::Function, propsd::Float64)
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
              propcov=nothing,
              propcov_init=eye(length(init)),
              temper_power=1,
              momentum=0.1, min_temper=1.0, max_temper=1.0,
              verbose=0) where T <: AbstractFloat
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
  if propcov == nothing
    propcov = PropCov(propcov_init, momentum=momentum, batchsize=batchsize)
  end

  # Update function
  function update(s, iter, propcov)
    # Compute temperature
    p = ((nburn - iter)^temper_power + 1) / (nburn^temper_power - 1)
    temper = iter < nburn ? max_temper^p : min_temper

    # Log density of target distribution
    _logprob(s::Vector{Float64}) = logprob(s, temper)

    # Update state
    s .= metropolis(s, _logprob, propcov)

    # Build a message
    msg = "temper: $(round(temper, digits=3))"

    return msg
  end


  # Run MCMC
  msg = ""
  for i in 1:total_iters
    print("\r$(i)/$(total_iters)")
    # Update state
    for t in 1:thin
      msg = update(state, i, propcov.value)
    end
    
    # Update proposal covariance
    update!(propcov, state)

    # Sanity check
    if verbose > 0 && mod(i, batchsize) == 0
      println(" | proposal covariance: $(round.(propcov.value[1], digits=3)) | $(msg)")
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
