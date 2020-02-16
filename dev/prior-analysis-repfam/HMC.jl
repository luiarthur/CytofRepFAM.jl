module HMC
using Distributions

"""
Get random momentum from a parameter.
"""
rand_momentum(x) = randn(size(x))


"""
Compute kinetic energy given momentum
"""
compute_kinetic_energy(momentum) = sum(momentum .^ 2) / 2

"""
See p. 14 of
https://arxiv.org/pdf/1206.1901.pdf
"""
# TODO: Check this
function hmc_update(curr_state,
                    logprob::Function, grad_logprob::Function,
                    num_leapfrog_steps::Int, eps::Float64;
                    momentum_sd::Real=1)
  # Position variables
  qs = copy(curr_state)

  # Momentum variables
  ps = rand_momentum(qs)

  # Current kinetic energy
  curr_K = compute_kinetic_energy(ps)
  
  # Potential enerdy
  U(s) = -logprob(s)

  # Current potential energy
  curr_U = U(qs)

  # Compute gradient of U
  grad_U(qs) = -grad_logprob(qs)

  # Make a half step for momentum at the beginning
  gs = grad_U(qs)
  ps .-= eps * qs / 2

  # Alternate full steps for position and momentum
  for i in 1:num_leapfrog_steps
    # Make a full step for the position
    qs .+= eps * ps

    # Make a full step for the momentum, except at end of trajectory
    if i < num_leapfrog_steps
      gs = grad_U(qs)
      ps .-= eps * qs
    end
  end

  # Make a half step for momentum at the end
  gs = grad_U(qs)
  ps .-= eps * qs / 2

  # Proposed potential energy
  cand_U = U(qs)

  # Proposed kinetic  energy
  cand_K = compute_kinetic_energy(ps)

  ### Metropolis step ###
  # log acceptance ratio.  NOTE: potential = -log_prob
  log_acceptance_ratio = curr_U + curr_K - cand_U - cand_K

  # Accept or reject
  if log_acceptance_ratio > log(rand())
    return qs, -cand_U
  else
    return curr_state, -curr_U
  end
  ### End of Metropolis step ###
end

end
