"""
Parallel Tempering
"""
module ParTemp

"""
Compute the log probability of swapping chains of different temperatures.

# Example

lpdf_normal(z) = -(log(2*pi) + z*z) / 2
loglike(s) = sum(lpdf_normal.(s))
states = (randn(100), randn(100))  # needs to be the same type
temps = (1.0, 2.0)  # each has to be greater than 1
compute_log_accept_ratio(loglike, states, temps)
compute_accept_ratio(loglike, states, temps)
"""
function compute_log_accept_ratio(loglike::Function, states, temps)
  s1, s2 = states
  t1, t2 = temps
  return min((1/t1 - 1/t2) * (loglike(s2) - loglike(s1)), 0)
end


"""
Compute the probability of swapping chains of different temperatures.

See doc for compute_log_accept_ratio
"""
function compute_accept_ratio(loglike::Function, states, temps)
  return exp(compute_log_accept_ratio(loglike, states, temps))
end

end  # ParTemp
