# TODO: Need to finish
#       Need to test
function fit_fs_tp!(init::StateFS,
                    cfs::ConstantsFS,
                    dfs::DataFS;
                    nmcmc::Int, nburn::Int, 
                    batchprop::Vector{Int}=0.1, prior_thin::Int=2,
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

  printMsg(iter::Int, msg::String) = if printFreq > 0 && iter % printFreq == 0
    print(msg)
  end

  @assert mod(num_tempers, 2) == 0
  if tfs != nothing
    @assert num_tempers == length(tfs)
  end

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
  function update(state, iter::Int)
    # Whether or not to marginalize over lambda and gamma.
    # We want to do this more often at the beginning, and less at the end.
    zmarg = ((Z_marg_lamgam - Z_marg_lamgam_min) * 
             exp(-iter/Z_marg_lamgam_decay_rate) + 
             Z_marg_lamgam_min) > rand()

    # Update state using trained prior
    update_via_trained_prior!(state, dfs, cfs, tfs, batchprop, prior_thin,
                              fix=fix, use_repulsive=use_repulsive,
                              Z_marg_lamgam=Z_marg_lamgam, sb_ibp=sb_ibp,
                              time_updates=time_updates)

    # Pull out inner componenets for convenience
    s = state.theta
    c = cfs.constants
    d = dfs.data

    # Append loglike
    append!(ll, compute_marg_loglike(s, c, d, 1.0))

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
  samples, state = MCMC.gibbs(deepcopy(init), update, monitors=monitors,
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
  out[:batchsize] = batchsize
  out[:prior_thin] = prior_thin

  return out
end
