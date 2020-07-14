function update_Z_v2!(s::State, c::Constants, d::Data, tuners::Tuners,
                      sb_ibp::Bool; use_repulsive::Bool=false)
  if 0.1 > rand()
    # update Z marginalizing over lam and gam
    update_Z_marg_lamgam!(s, c, d, sb_ibp, use_repulsive=use_repulsive)
  else
    update_Z!(s, c, d, sb_ibp, use_repulsive=use_repulsive)
  end
end

function update_Z_marg_lamgam!(j::Int, k::Int,
                               A::Vector{Vector{Float64}},
                               B0::Vector{Matrix{Float64}},
                               B1::Vector{Matrix{Float64}},
                               s::State, c::Constants, d::Data, sb_ibp::Bool;
                               use_repulsive::Bool=false)
  v = sb_ibp ? cumprod(s.v) : s.v
  Z0 = deepcopy(s.Z)
  Z0[j, k] = false 
  lp0 = MCMC.log1m(v[k]) + log_dmix_nolamgam(Z0, A, B0, B1, s, c, d)

  Z1 = deepcopy(s.Z)
  Z1[j, k] = true 
  lp1 = log(v[k]) + log_dmix_nolamgam(Z1, A, B0, B1, s, c, d)

  if use_repulsive
    lp0 += log_penalty_repFAM(k, Z0, c.log_repulsive_fn)
    lp1 += log_penalty_repFAM(k, Z1, c.log_repulsive_fn)
  end

  p1_post = 1 / (1 + exp(lp0 - lp1))
  new_Zjk_is_one = p1_post > rand()
  s.Z[j, k] = new_Zjk_is_one
end


function update_Z_marg_lamgam!(s::State, c::Constants, d::Data, sb_ibp::Bool;
                               use_repulsive::Bool=false)
  # Precompute mus0, mus1, sig
  mus0 = mus(false, s, c, d)
  mus1 = mus(true, s, c, d)
  sig = sqrt.(s.sig2 * c.temper)

  # Precompute A, B0, B1

  # size(A[i]) = N[i]
  A = [if s.eps[i] > 0
         [logdnoisy(i, n, s, c, d) for n in 1:d.N[i]]
       else
         zeros(d.N[i])
       end for i in 1:d.I]

  # size(B0[i]) = (N[i], J)
  B0 = [[logdmixture(s.y_imputed[i][n, j], mus0, sig[i], s.eta[false][i, j, :])
         for n in 1:d.N[i], j in 1:d.J] for i in 1:d.I]

  # size(B1[i]) = (N[i], J)
  B1 = [[logdmixture(s.y_imputed[i][n, j], mus1, sig[i], s.eta[true][i, j, :])
         for n in 1:d.N[i], j in 1:d.J]
        for i in 1:d.I]

  # NOTE: The updates of Z_{j,k} cannot be made asynchronous naively.
  for j in 1:d.J
    for k in 1:c.K
      update_Z_marg_lamgam!(j, k, A, B0, B1, s, c, d, sb_ibp,
                            use_repulsive=use_repulsive)
    end
  end
end
