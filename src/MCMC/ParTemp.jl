"""
Parallel Tempering
"""
module ParTemp

function compute_log_accept_ratio(loglike::Function, states, temps)
  s1, s2 = states
  t1, t2 = temps
  return (1/t1 - 1/t2) * (loglike(s2) - loglike(s1))
end

function compute_accept_ratio(loglike::Function, states, temps)
  return exp(compute_log_accept_ratio(loglike, states, temps))
end

end  # ParTemp
