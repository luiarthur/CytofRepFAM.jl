function update_omega!(s::StateFS, c::ConstantsFS, d::DataFS, t::TunersFS)
  Q = length(s.omega)
  for q in 1:Q
    update_omega!(q, s, c, d, t)
  end
end

function update_omega!(q::Integer, s::StateFS, c::ConstantsFS, d::DataFS,
                       t::TunersFS)

  function loglike(om_q::Float64)::Float64
    om = copy(s.omega)  # make a copy
    om[q] = om_q
    p = compute_p(d.X, om)

    ll = 0.
    for i in 1:d.data.I
      lpdf(rik::Bool)::Float64 = logpdf(Bernoulli(p[i]), rik)
      for k in 1:c.constants.K
        ll += lpdf(s.r[i, k])
      end
    end

    return ll
  end

  logprior(om_q::Float64)::Float64 = logpdf(c.omega_prior, om_q)
  log_prob(om_q::Float64)::Float64 = loglike(om_q) + logprior(om_q)

  s.omega[q] = MCMC.metropolisAdaptive(s.omega[q], log_prob, t.omega[q])
end

# Omega = p when using Cell means model for X.
function update_omega_cellmeans!(s::StateFS, c::ConstantsFS, d::DataFS,
                                 t::TunersFS)
  Q = length(s.omega)
  for q in 1:Q
    update_omega_cellmeans!(q, s, c, d, t)
  end
end

function update_omega_cellmeans!(q::Integer, s::StateFS, c::ConstantsFS,
                                 d::DataFS, t::TunersFS)
  a_new, b_new = params(c.p_prior) 

  xi = findfirst(z -> z == 1, d.X[i, :])
  K = c.constants.K
  for i in 1:d.data.I
    if xi == q
      ri_sum = sum(s.r[i, :])
      a_new += ri_sum
      b_new += (K - ri_sum)
    end
  end

  s.omega[q] = rand(Beta(a_new, b_new))
end
