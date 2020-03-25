# TODO: Need to finish
#       Need to test
function fit_fs_tp!(init::StateFS,
                    cfs::ConstantsFS,
                    dfs::DataFS;
                    nmcmc::Int, nburn::Int, 
                    # These args are for iMCMC:
                    batchprop::Float64=0.1, prior_thin::Int=2,
                    temper::Float64=1.0, mb_update_burn_prop=0.6, anneal=false,
                    # End of iMCMC args.
                    tfs::Union{Nothing, Vector{TunersFS}}=nothing,
                    monitors=[monitor1, monitor2],
                    fix::Vector{Symbol}=Symbol[],
                    thins::Vector{Int}=[2, nsamps_to_thin(10, nmcmc)],
                    thin_dden::Int=1,
                    printFreq::Int=0, 
                    checkpoint=0,
                    computeDIC::Bool=false, computeLPML::Bool=false,
                    computedden::Bool=false,
                    sb_ibp::Bool=false,
                    use_repulsive::Bool=true,
                    Z_marg_lamgam::Float64=1.0,
                    Z_marg_lamgam_decay_rate::Float64=100.0,
                    Z_marg_lamgam_min::Float64=0.05,
                    verbose::Int=1,
                    time_updates=false,
                    seed::Int=-1)

  @assert 0 <= mb_update_burn_prop <= 1

  printMsg(iter::Int, msg::String) = if printFreq > 0 && iter % printFreq == 0
    print(msg)
  end

  # Set random seed if needed
  if seed >= 0
    Random.seed!(seed)
  end

  # Print temperature
  println("temperature: $(temper)")

  # Create Tuners if needed 
  if tfs == nothing
    tfs = TunersFS(Tuners(dfs.data.y, cfs.constants.K),
                   init.theta,
                   dfs.X)
  end

  if verbose >= 1
    fixed_vars_str = join(fix, ", ")
    if fixed_vars_str == ""
      fixed_vars_str = "nothing"
    end
    println("fixing: $fixed_vars_str")
    println("Use stick-breaking IBP: $(sb_ibp)")
    println("Z_marg_lamgam: $(Z_marg_lamgam)")
    println("Z_marg_lamgam_decay_rate: $(Z_marg_lamgam_decay_rate)")
    println("Z_marg_lamgam_min: $(Z_marg_lamgam_min)")
    println("use_repulsive: $(use_repulsive)")
    flush(stdout)
  end

  @assert printFreq >= -1
  if printFreq == 0
    numPrints = 10
    printFreq = Int(ceil((nburn + nmcmc) / numPrints))
  end

  # Instantiate (but not initialize) CPO stream
  if computeLPML
    cpoStream = MCMC.CPOstream{Float64}()
  end

  # DIC
  if computeDIC
    dicStream = initDicStream(dfs.data)
    loglikeDIC(param::DICparam) = computeLoglikeDIC(dfs.data, param)

    convertStateToDicParam(s::State)::DICparam = let
      _convertStateToDicParam(s, cfs.constants, dfs.data)
    end
  end

  # Cache for data density
  dden = Matrix{Vector{Float64}}[]

  # Cache for loglikelihood
  ll = Float64[]

  # Update function
  function update!(state, iter::Int, out)
    # Whether or not to marginalize over lambda and gamma.
    # We want to do this more often at the beginning, and less at the end.
    zmarg = ((Z_marg_lamgam - Z_marg_lamgam_min) * 
             exp(-iter/Z_marg_lamgam_decay_rate) + 
             Z_marg_lamgam_min) > rand()

    # Iterations to do minibatch update
    iters_mb_update_all_params = round(Int, nburn * mb_update_burn_prop)

    # Calculate current temperature
    curr_temper = if anneal
      # Gradually cool the temperature so that it is 1 when
      # we reach iteration `iters_mb_update_all_params`, which
      # is less than `nburn`.
      itmax = iters_mb_update_all_params
      power = (itmax - iter + 1) / itmax
      clamp(temper^power, 1.0, Inf)
    else
      temper
    end
    if verbose > 0
      print(" -- Current temper: $(curr_temper)")
    end

    # Update state using trained prior
    minibatch_update_all_params = iter < iters_mb_update_all_params
    update_via_trained_prior!(state, dfs, cfs, tfs, batchprop, prior_thin,
                              fix=fix, use_repulsive=use_repulsive,
                              Z_marg_lamgam=zmarg, sb_ibp=sb_ibp,
                              time_updates=time_updates, temper=curr_temper,
                              minibatch_update_all_params=minibatch_update_all_params,
                              verbose=verbose-2)

    # Pull out inner componenets for convenience
    s = state.theta
    c = cfs.constants
    d = dfs.data

    # Append loglike
    append!(ll, compute_marg_loglike(s, c, d, 1.0) / curr_temper)

    if computedden && iter > nburn && (iter - nburn) % thin_dden == 0
      # NOTE: `datadensity(s, c, d)` returns an (I x J) matrix of vectors of
      # length g.
      append!(dden, [datadensity(s, c, d)])
    end

    if computeLPML && iter > nburn
      # Inverse likelihood for each data point
      like = [[compute_like(i, n, s, c, d)
               for n in 1:d.N[i]] for i in 1:d.I]

      # Update (or initialize) CPO
      MCMC.updateCPO(cpoStream, vcat(like...))

      # Add to printMsg
      printMsg(iter, " -- LPML: $(MCMC.computeLPML(cpoStream))")
    end

    if computeDIC && iter > nburn
      # Update DIC
      MCMC.updateDIC(dicStream, s, updateParams,
                     loglikeDIC, convertStateToDicParam)

      # Add to printMsg
      printMsg(iter, " -- DIC: $(MCMC.computeDIC(dicStream, loglikeDIC,
                                                 paramMeanCompute))")

      DICg = MCMC.DIC_gelman(ll[(nburn+1):end])
      printMsg(iter, " -- DIC_gelman: $(DICg)")
    end

    printMsg(iter, "\n")
    flush(stdout)
  end

  println("Running Gibbs sampler ...")
  samples, state = MCMC.gibbs(deepcopy(init), update!, monitors=monitors,
                              thins=thins, nmcmc=nmcmc, nburn=nburn,
                              printFreq=printFreq, loglike=ll,
                              printlnAfterMsg=false)

  out = Dict(:samples => samples,
             :lastState => state,
             :init => init,
             :ll => ll,
             :c => cfs,
             :d => dfs)

  if computeDIC || computeLPML
    LPML = computeLPML ? MCMC.computeLPML(cpoStream) : NaN
    Dmean, pD = computeDIC ? MCMC.computeDIC(dicStream, loglikeDIC,
                                             paramMeanCompute,
                                             return_Dmean_pD=true) : (NaN, NaN)

    DICg = MCMC.DIC_gelman(ll[(nburn+1):end])

    metrics = Dict(:LPML => LPML,
                   :DIC => Dmean + pD,
                   :Dmean => Dmean,
                   :pD => pD,
                   :DICg => DICg)
    println()
    println("metrics:")
    for (k, v) in metrics
      out[k] = v
      println("$k => $v")
    end
    flush(stdout)
  end

  if computedden
    out[:dden] = dden
  end

  out[:nburn] = nburn
  out[:nmcmc] = nmcmc
  out[:batchprop] = batchprop
  out[:prior_thin] = prior_thin
  out[:temper] = temper
  out[:anneal] = anneal

  return out
end
