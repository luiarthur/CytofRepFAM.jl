mutable struct Metropolis{T <: AbstractFloat, S}
  state::S
  logprob_state::T
end

function update!(m::Metropolis{T, S},
                 logprob::Function,
                 logjump::Function,
                 randjump::Function;
                 verbose::Int=0) where {T <: AbstractFloat, S}
  # Draw candidate from jump (given current state)
  cand = randjump(m.state)

  # Compute log density of candidate draw 
  logprob_cand = logprob(cand)

  # Compute log porposal density of prev -> candidate
  logjump_cand_from_prev = logjump(cand, m.state)

  # Compute log porposal density of candidate -> prev
  logjump_prev_from_cand = logjump(m.state, cand)

  # Compute log acceptance ratio
  log_accept_ratio = ((logprob_cand - m.logprob_state) -
                      (logjump_cand_from_prev - logjump_prev_from_cand))

  # Log of rand uniform
  log_u = log(rand())

  # With some probability, update current state with proposed state.
  if log_accept_ratio > log_u
    if verbose > 0
      alpha = round(log_accept_ratio, digits=3)
      println("Accepted proposal with log_accept_ratio = $(alpha)")
    end
    m.state = deepcopy(cand)
    m.logprob_state = logprob_cand
  elseif verbose > 0
    alpha = round(log_accept_ratio, digits=3)
    println("Rejected proposal with log_accept_ratio = $(alpha)")
  end

  return m
end
