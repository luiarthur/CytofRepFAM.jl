function update_r!(s::StateFS, c::ConstantsFS, d::DataFS)
  # update each r_{i, k}
  for i in 1:d.data.I
    ldmix = logdmixture_r(i, s.theta, c.constants, d.data)
    for k in 1:c.constants.K
      update_r_marg_lam!(i, k, ldmix, s, c, d)
      # NOTE: Make sure to always update W immediately after updating r or
      #       W_star!
      update_W!(s, c, d)
    end
  end
end


function compute_w(ws::Vector{Float64}, r::Vector{Bool})
  wr = ws .* r
  return wr ./ sum(wr)
end


function logdmixture_r(i::Integer, s::State, c::Constants, d::Data)
  """
  Precomputes the relevant mixture density, marginalized over gam & lam.
  This is possible because (Z, mus, sig, eta, y_miss) are held constant.
  Note that W is a function of (r, W*), and hence is NOT held constant, 
  and so does not appear here. 
  """
  # Precompute / reshapes
  mus0 = reshape(mus(false, s, c, d), 1, 1, c.L[0])
  mus1 = reshape(mus(true, s, c, d), 1, 1, c.L[1])
  sig = sqrt.(s.sig2 * c.temper)
  Ni = d.N[i]
  Z = reshape(s.Z, 1, d.J, c.K)
  eta0i = s.eta[false][i:i, :, :]
  eta1i = s.eta[true][i:i, :, :]

  # Ni x J x 1
  yi = reshape(s.y_imputed[i], Ni, d.J, 1)

  # Ni x J x 1
  kernel0 = MCMC.lpdf_gmm(yi, mus0, sig[i], eta0i, dims=3)
  kernel1 = MCMC.lpdf_gmm(yi, mus1, sig[i], eta1i, dims=3)

  # (1 x J x K) .* (Ni x 1 x K) -> Ni x K
  Z_mix = let
    tmp = sum(Z .* kernel1 + (1 .- Z) .* kernel0, dims=2)
    dropdims(tmp, dims=2)
  end

  return Z_mix
end


function update_r_marg_lam!(i::Integer, k::Integer, ldmix::Matrix{Float64},
                            s::StateFS, c::ConstantsFS, d::DataFS)

  # Log likelihood as a function of `r_{i, k}`
  function loglike(r_ik::Bool)::Float64
    # Make a copy of the current (k-dim) vector r_i
    ri = deepcopy(s.r[i, :])

    # Replace r_{i, k} with the provided `r_ik`
    ri[k] = r_ik

    # Compute W given W* and the updated r_i
    wi = compute_w(s.W_star[i, :], ri)

    # Get the indices of the W_i such that W_{i, k} > 0
    non_zero_wri_idx = findall(ri)

    # Initialize log likelihood
    ll = 0.0

    if length(non_zero_wri_idx) > 0  # i.e. if any W_i > 0
      wi_rs = reshape(wi[non_zero_wri_idx], 1, sum(ri))
      lln = MCMC.logsumexp(log.(wi_rs) .+ ldmix[:, non_zero_wri_idx], dims=2)
      ll = sum(lln)
    else
      ll = -Inf
    end

    return ll
  end

  # Log prior as a function of `r_{i, k}`
  function logprior(r_ik::Bool)::Float64
    p_xi = compute_p(d.X[i, :], s.omega)
    return logpdf(Bernoulli(p_xi), r_ik)
  end

  # Log full conditional as a function of `r_{i, k}`
  log_prob(r_ik::Bool)::Float64 = logprior(r_ik) + loglike(r_ik)

  # update r_{i, k} with a metropolis step
  rik_metropolis_update!(i, k, log_prob, s)
end


function rik_metropolis_update!(i::Integer, k::Integer,
                                log_prob::Function, s::StateFS)
  cand = !s.r[i, k]  # Flip bit
  log_acceptance_ratio = log_prob(cand) - log_prob(s.r[i, k])

  # accept with probability `accecptance_ratio`
  if log_acceptance_ratio > log(rand())
    s.r[i, k] = cand
  end
end
