# FIXME! 
# Need to update.
function fit_fs_pt!(cfs::ConstantsFS,
                    dfs::DataFS;
                    tempers::Vector{Float64},
                    nmcmc::Int, nburn::Int, 
                    tfs::Union{Nothing, Vector{TunersFS}}=nothing,
                    inits=nothing,
                    return_all_states=false,
                    remove_current_workers::Bool=true,
                    monitors=[monitor1, monitor2],
                    fix::Vector{Symbol}=Symbol[],
                    thins::Vector{Int}=[2, nsamps_to_thin(10, nmcmc)],
                    thin_dden::Int=1,
                    printFreq::Int=0, 
                    randpair=0.0,
                    checkpoint=0,
                    computeDIC::Bool=false, computeLPML::Bool=false,
                    computedden::Bool=false,
                    sb_ibp::Bool=false,
                    use_repulsive::Bool=true,
                    Z_marg_lamgam::Float64=1.0,
                    Z_marg_lamgam_decay_rate::Float64=100.0,
                    Z_marg_lamgam_min::Float64=0.05,
                    verbose::Int=1,
                    swap_freq::Float64=1.0,
                    time_updates=false,
                    seed::Int=-1)

  printMsg(iter::Int, msg::String) = if printFreq > 0 && iter % printFreq == 0
    print(msg)
  end

  # Number of temperatures
  num_tempers = length(tempers)

  @assert mod(num_tempers, 2) == 0
  if tfs != nothing
    @assert num_tempers == length(tfs)
  end

  # Assert swap frequency is in unit interval
  @assert 0 < swap_freq <= 1.0

  # Cache for swap counts
  swapcounts = zeros(Int, num_tempers, num_tempers)

  # Cache for pair counts
  paircounts = zeros(Int, num_tempers, num_tempers)

  # Set random seed if needed
  if seed >= 0
    Random.seed!(seed)
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

  ### Coldest chain stuff ###

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

  # Update function
  # function update(states::Vector{StateFS}, args::Vector{Dict{Symbol, Any}}, iter::Int)
  function update(states, args, iter::Int)
    # Whether or not to marginalize over lambda and gamma.
    # We want to do this more often at the beginning, and less at the end.
    zmarg = ((Z_marg_lamgam - Z_marg_lamgam_min) * 
             exp(-iter/Z_marg_lamgam_decay_rate) + 
             Z_marg_lamgam_min) > rand()

    s_arg_vec = pmap(states, args) do s, arg
      tu = (arg[:c].constants.temper == 1) && time_updates
      
      # For reproducibility.
      if seed > - 1
        temper_idx = findfirst(tau -> tau == arg[:c].constants.temper, tempers)
        Random.seed!(iter + temper_idx + seed)
      end

      update_state_feature_select!(s, arg[:c], dfs, arg[:t],
                                   ll=arg[:ll], fix=fix,
                                   use_repulsive=use_repulsive,
                                   # TODO: make this more explicit, instead
                                   # of random. See `update_Z_v2`.
                                   Z_marg_lamgam=zmarg,
                                   sb_ibp=sb_ibp, time_updates=tu)
      (s, arg)
    end

    states = [sa[1] for sa in s_arg_vec]
    args = [sa[2] for sa in s_arg_vec]

    # Swap Chains
    if verbose > 0
      println()
    end
    llf(s, tuner) = compute_marg_loglike(s, cfs.constants, dfs.data, tuner)

    if swap_freq > rand()
      swapchains!(states, llf, tempers,
                  paircounts=paircounts, swapcounts=swapcounts,
                  randpair=randpair, verbose=verbose)
    end
    if iter == nburn
      println("swapcounts / paircounts:")
      println(swapcounts ./ (paircounts .+ 1e-6))
      println("Resetting swapcounts, paircounts ...")
      swapcounts .= 0.0
      paircounts .= 0.0
    end

    s = states[1]
    c = args[1][:c]
    d = dfs
    ll = args[1][:ll]

    if computedden && iter > nburn && (iter - nburn) % thin_dden == 0
      # NOTE: `datadensity(s, c, d)` returns an (I x J) matrix of vectors of
      # length g.
      append!(dden, [datadensity(s.theta, c.constants, d.data)])
    end

    if computeLPML && iter > nburn
      # Inverse likelihood for each data point
      like = [[compute_like(i, n, s.theta, c.constants, d.data)
               for n in 1:d.data.N[i]] for i in 1:d.data.I]

      # Update (or initialize) CPO
      MCMC.updateCPO(cpoStream, vcat(like...))

      # Add to printMsg
      printMsg(iter, " -- LPML: $(MCMC.computeLPML(cpoStream))")
    end

    if computeDIC && iter > nburn
      # Update DIC
      MCMC.updateDIC(dicStream, s.theta, updateParams,
                     loglikeDIC, convertStateToDicParam)

      # Add to printMsg
      printMsg(iter, " -- DIC: $(MCMC.computeDIC(dicStream, loglikeDIC,
                                                 paramMeanCompute))")

      DICg = MCMC.DIC_gelman(ll[(nburn+1):end])
      printMsg(iter, " -- DIC_gelman: $(DICg)")
    end

    printMsg(iter, "\n")
    flush(stdout)

    return states, args
  end


  # Create vectors of states
  states = if inits == nothing
    println("Using random inits!")
    [let
       s = genInitialState(cfs.constants, dfs.data)
       s.eps .= 0.0
       StateFS{Float64}(s, dfs)
     end for _ in tempers]
  else
    @assert length(inits) == num_tempers
    inits
  end

  # Create Args
  args = [let
            ll = Float64[]
            c = deepcopy(cfs)
            c.constants.temper = tempers[i]
            t = if tfs == nothing
              TunersFS(Tuners(dfs.data.y, cfs.constants.K),
                       states[1].theta, dfs.X)
            else
              tfs[i]
            end
            Dict(:ll => ll, :c => c, :t => t)
          end
          for i in 1:num_tempers]

  println("Running Gibbs sampler ...")
  samples, states, args = MCMC.gibbs_pt(states, args, update, monitors=monitors,
                                        thins=thins, nmcmc=nmcmc, nburn=nburn,
                                        printFreq=printFreq,
                                        printlnAfterMsg=false)

  out = Dict(:samples => samples,
             :lastState => states[1],
             :lls => [arg[:ll] for arg in args],
             :c => cfs,
             :d => dfs,
             :tempers => tempers,
             :swap_freq => swap_freq,
             :randpair => randpair,
             :paircounts => paircounts,
             :swapcounts => swapcounts)

  if computeDIC || computeLPML
    LPML = computeLPML ? MCMC.computeLPML(cpoStream) : NaN
    Dmean, pD = computeDIC ? MCMC.computeDIC(dicStream, loglikeDIC,
                                             paramMeanCompute,
                                             return_Dmean_pD=true) : (NaN, NaN)

    ll1 = args[1][:ll]
    DICg = MCMC.DIC_gelman(ll1[(nburn+1):end])

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

  # Return all states if requested
  if return_all_states
    out[:all_last_states] = deepcopy(states)
  end

  return out
end
