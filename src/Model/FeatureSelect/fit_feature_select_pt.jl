function swap!(i, j, x)
  tmp = deepcopy(x[i])
  x[i] = deepcopy(x[j])
  x[j] = tmp
  return x
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

    s1, s2 = states[i].theta, states[j].theta
    t1, t2 = temperatures[i], temperatures[j]
    log_accept_ratio = MCMC.WSPT.compute_log_accept_ratio(loglike,
                                                          (s1, s2),
                                                          (t1, t2))
    should_swap_chains = log_accept_ratio > log(rand())

    if verbose > 2
      println("    mu*0 for state $(i): $(-cumsum(states[i].theta.delta[false]))")
      println("    mu*1 for state $(i): $(cumsum(states[i].theta.delta[true]))")
      println("    sig2 for state $(i): $(states[i].theta.sig2)")
      println("    mu*0 for state $(j): $(-cumsum(states[j].theta.delta[false]))")
      println("    mu*1 for state $(j): $(cumsum(states[j].theta.delta[true]))")
      println("    sig2 for state $(j): $(states[j].theta.sig2)")
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

function fit_fs_pt!(init::StateFS, cfs::ConstantsFS, dfs::DataFS, tfs::TunersFS;
                    use_rand_inits=false,
                    tempers::Vector{Float64}, ncores::Int,
                    nmcmc::Int, nburn::Int, 
                    swap_freq::Int=10,
                    remove_current_workers::Bool=true,
                    monitors=[monitor1, monitor2],
                    fix::Vector{Symbol}=Symbol[],
                    thins::Vector{Int}=[2, nsamps_to_thin(10, nmcmc)],
                    thin_dden::Int=1,
                    printFreq::Int=0, 
                    computeDIC::Bool=false, computeLPML::Bool=false,
                    computedden::Bool=false,
                    sb_ibp::Bool=false,
                    use_repulsive::Bool=true, Z_marg_lamgam::Bool=true,
                    verbose::Int=1,
                    time_updates=[false for _ in tempers],
                    Z_thin::Int=1, seed::Int=-1)

  # Number of temperatures
  num_tempers = length(tempers)
  @assert mod(num_tempers, 2) == 0

  # Create distributed arrays for:
  # - loglike, ConstantsFS, StatesFS, TunersFS
  lls = [Float64[0.0] for _ in tempers]
  cs = [let 
          ct = deepcopy(cfs)
          ct.constants.temper = tau
          ct
        end for tau in tempers]
  states = if use_rand_inits
    println("Using random inits!")
    [let
       s = genInitialState(cfs.constants, dfs.data)
       StateFS{Float64}(s, dfs)
     end for _ in tempers]
  else
    [deepcopy(init) for _ in tempers]
  end
  tuners = [deepcopy(tfs) for _ in tempers]
  swapcounts = zeros(Int, num_tempers, num_tempers)
  paircounts = zeros(Int, num_tempers, num_tempers)

  println("About to start parallel chains ...")
  for iter in 1:(nburn+nmcmc)
    # Perform updates in parallel
    # @time out = pmap(states, cs, tuners, lls, time_updates) do s, c, t, ll, tu
    #   @time for _ in 1:swap_freq
    out = pmap(states, cs, tuners, lls, time_updates) do s, c, t, ll, tu
      for _ in 1:swap_freq
        update_state_feature_select!(s, c, dfs, t,
                                     ll=ll, fix=fix,
                                     use_repulsive=use_repulsive,
                                     # TODO: make this more explicit, instead
                                     # of random. See `update_Z_v2`.
                                     Z_marg_lamgam=Z_marg_lamgam,
                                     randpair=0.0,
                                     sb_ibp=sb_ibp, time_updates=tu)
        # println(ll)
      end
      return Dict(:s => s, :c => c, :t => t, :ll => ll)
    end
    println("Done with iter $iter")

    # replace stuff
    states = [o[:s] for o in out]
    cs = [o[:c] for o in out]
    tuners = [o[:t] for o in out]
    lls = [o[:ll] for o in out]

    if iter > nburn  # && iter % swap_freq == 0
      llf(s, t) = compute_marg_loglike(s, cfs.constants, dfs.data, t)
      swapchains!(states, llf, tempers,
                  paircounts=paircounts, swapcounts=swapcounts,
                  randpair=randpair,
                  verbose=verbose)
    end
  end
  println("After running parallel chains ...")

  out = Dict(:lls => lls,
             :state => states[1],
             :c => cfs,
             :d => dfs,
             :t => tfs,
             :paircounts => paircounts,
             :swapcounts => swapcounts)
  return out
end


#= Possible implementation detail
@assert mod(nmcmc + nburn, swap_freq) == 0

for iter in 1:div(nmcmc + nburn, swap_freq)
  # Update states in parallel
  for _ in 1:swap_freq
    states = update(states)
  end

  # SWAP STATES
end
=#
