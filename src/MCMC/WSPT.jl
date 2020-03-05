"""
Weight Stabilizing Parallel Tempering

See: https://link.springer.com/article/10.1007/s11222-019-09863-3
"""
module WSPT

function gentempers(maxtemp, ntemps; degree=1)
  i = collect(1:ntemps)
  powers = (i.^ degree .- 1) / (ntemps ^ degree - 1)
  return maxtemp .^ powers
end

function compute_log_accept_ratio(loglike::Function, states, temps; verbose=0)
  s1, s2 = states
  t1, t2 = temps

  ll_s1_t2 = loglike(s1, t2)
  ll_s2_t1 = loglike(s2, t1)
  ll_s1_t1 = loglike(s1, t1)
  ll_s2_t2 = loglike(s2, t2)

  if verbose > 0
    println("ll_s1_t2: $(ll_s1_t2)")
    println("ll_s2_t1: $(ll_s2_t1)")
    println("ll_s1_t1: $(ll_s1_t1)")
    println("ll_s2_t2: $(ll_s2_t2)")
  end

  loglike_cand = ll_s1_t2 + ll_s2_t1
  loglike_curr = ll_s1_t1 + ll_s2_t2

  return loglike_cand - loglike_curr
end

function compute_accept_ratio(loglike::Function, states, temps; verbose=0)
  return exp(compute_log_accept_ratio(loglike, states, temps, verbose=verbose))
end

end
