# TODO: Implement the log-likelihood, which marginalizes over lam, gam.
# This will be used for parallel tempering.
"""
Basing this off of conditional DIC:
p(y_in, m_in | theta, beta) \\propto
    p(y_obs | y_mis, theta) * p(m | y_obs, y_mis, beta)
We don't have this exactly when we marginalize out the indicators,
but we can have

p(y_complete, m | theta)

where y_complete is the (y_obs, y_mis), for a given iteration.
"""
function compute_marg_loglike(s::State, c::Constants, d::Data, temper::Float64;
                              y=nothing, compute_prob_miss=true)
  # Boolean: Whether y was supplied
  y_supplied = (y != nothing)

  if !y_supplied
    y = s.y_imputed
  end

  # log-likelihood
  ll = 0.0

  # Dimensions
  J = d.J
  K = c.K
  L = c.L

  # Precompute / reshapes
  mus0 = reshape(mus(false, s, c, d), 1, 1, L[0])
  mus1 = reshape(mus(true, s, c, d), 1, 1, L[1])
  sig = sqrt.(s.sig2)

  for i in 1:d.I
    # Dimensions
    Ni = if y_supplied
      size(y[i], 1)
    else
      d.N[i]
    end
    Z = reshape(s.Z, 1, J, K)
    eta0i = s.eta[false][i:i, :, :]
    eta1i = s.eta[true][i:i, :, :]

    # Ni x J x 1
    yi = reshape(y[i], Ni, J, 1)

    # Ni x J x 1
    kernel0 = MCMC.lpdf_gmm(yi, mus0, sig[i], eta0i, temper, dims=3)
    kernel1 = MCMC.lpdf_gmm(yi, mus1, sig[i], eta1i, temper, dims=3)

    # (1 x J x K) .* (Ni x 1 x K) -> Ni x K
    Z_mix = let
      tmp = sum(Z .* kernel1 + (1 .- Z) .* kernel0, dims=2)
      dropdims(tmp, dims=2)
    end

    # (Ni x K) .+ (1 x K) --> (Ni x K)
    f = Z_mix .+ log.(s.W[i:i, :])

    # Ni-dim
    lli = MCMC.logsumexp(f, dims=2)

    if compute_prob_miss
      # Add probability of missing for only imputed values only, because for
      # observed values, the probability evaluates to a constant across models.
      mi = Bool.(d.m[i])
      pm = [Model.prob_miss(y_inj, c.beta[:, i]) for y_inj in s.y_imputed[i][mi]]

      # Increment loglikelihood with missing prob stuff
      ll += sum(log.(pm))
    end

    # Increment loglikelihood
    ll += sum(lli)
  end  # i loop

  return ll
end
