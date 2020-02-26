"""
Weight Stabilizing Parallel Tempering

See: https://link.springer.com/article/10.1007/s11222-019-09863-3
"""
module WSPT

function compute_log_accept_ratio(loglike::Function, states, temps)
  s1, s2 = states
  t1, t2 = temps

  loglike_cand = loglike(s1, t2) + loglike(s2, t1)
  loglike_curr = loglike(s1, t1) + loglike(s2, t2)

  return min(loglike_cand - loglike_curr, 0)
end

function compute_accept_ratio(loglike::Function, states, temps)
  return exp(compute_log_accept_ratio(loglike, states, temps))
end

end
