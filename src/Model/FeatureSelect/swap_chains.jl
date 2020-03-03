function swap!(i, j, x)
  tmp = deepcopy(x[i])
  x[i] = deepcopy(x[j])
  x[j] = tmp
  return x
end

function println_mu_sig2(states, i)
  println("    mu*0 for state $(i): $(-cumsum(states[i].theta.delta[false]))")
  println("    mu*1 for state $(i): $(cumsum(states[i].theta.delta[true]))")
  println("    sig2 for state $(i): $(states[i].theta.sig2)")
end

function swapchains!(states, loglike, temperatures;
                     paircounts, swapcounts, verbose=0, randpair=0.0)
  """
  randpair: The proportion of time to propose swapping random pairs. 
            If set to 0.0, then pairs are swapped from hottest to coolest
            or coolest to hottest. If set to 1.0, then random pairs are 
            formed and proposed to be swapped.

  See: https://academic.oup.com/gji/article/196/1/357/585739
  """
  nchains = length(states)
  @assert (nchains == length(temperatures)) && (mod(nchains, 2) == 0)

  pairs = if randpair > rand()
    Iterators.partition(Random.shuffle(1:nchains), 2)
  else
    p = zip(1:(nchains - 1), 2:nchains)
    rand(Bool) ? p : Iterators.reverse(p)
  end

  for (i, j) in pairs
    # Increment pair counts
    paircounts[i, j] += 1
    paircounts[j, i] += 1

    si, sj = states[i].theta, states[j].theta
    ti, tj = temperatures[i], temperatures[j]
    log_accept_ratio = MCMC.WSPT.compute_log_accept_ratio(loglike,
                                                          (si, sj),
                                                          (ti, tj))
    should_swap_chains = log_accept_ratio > log(rand())

    if verbose > 2
      println_mu_sig2(states, i)
      println_mu_sig2(states, j)
    elseif verbose > 1
      println("swap log accept ratio ($(i), $(j)): $(log_accept_ratio)")
    end

    if should_swap_chains
      swap!(i, j, states)
      # Increment swap-counts matrix
      swapcounts[i, j] += 1
      swapcounts[j, i] += 1
      if verbose > 0
        println("Swapped chains $(i) and $(j)")
      end
    end
  end

  return
end
