nsamps_to_thin(nsamps::Int, nmcmc::Int) = max(1, div(nmcmc, nsamps))

monitor1 = [:theta__Z, :theta__v, :theta__alpha,
            :omega, :r, :theta__lam, :W_star, :theta__eta,
            :theta__W, :theta__delta, :theta__sig2]

monitor2 = [:theta__y_imputed, :theta__gam]


"""
printFreq: defaults to 0 => prints every 10%. turn off printing by 
           setting to -1.
thins: defaults to `[2, nsamps_to_thin(10, nmcmc)]`.
       i.e. Keep every other sample for monitor 1 (so the total number of
       samples is `nmcmc/2`). And keep 10 samples total for monitor 2.
"""
function fit_fs!(init::StateFS, c::ConstantsFS, d::DataFS;
                 nmcmc::Int, nburn::Int, 
                 tuners::Union{Nothing, TunersFS}=nothing,
                 monitors=[monitor1, monitor2],
                 fix::Vector{Symbol}=Symbol[],
                 thins::Vector{Int}=[2, nsamps_to_thin(10, nmcmc)],
                 thin_dden::Int=1,
                 printFreq::Int=0, 
                 computeDIC::Bool=false, computeLPML::Bool=false,
                 computedden::Bool=false,
                 sb_ibp::Bool=false,
                 use_repulsive::Bool=true, Z_marg_lamgam::Bool=true,
                 verbose::Int=1, time_updates::Bool=false, Z_thin::Int=1,
                 seed::Int=-1)

  # Set random seed if needed
  if seed >= 0
    Random.seed!(seed)
  end

  # We don't want to use noisy distribution.
  # Assert that all eps == 0
  if any(init.theta.eps .!= 0)
    msg = "WARNING: `init.theta.eps`: $(init.theta.eps) contains non-zeros values! "
    msg *= "Setting all eps to 0!"
    println(msg)

    # Set eps == 0
    init.theta.eps .= 0
  end
  flush(stdout)

  @assert Z_thin >= 0
  if Z_thin == 0
    Z_thin = d.data.J
  end

  if verbose >= 1
    fixed_vars_str = join(fix, ", ")
    if fixed_vars_str == ""
      fixed_vars_str = "nothing"
    end
    println("fixing: $fixed_vars_str")
    println("Use stick-breaking IBP: $(sb_ibp)")
    println("Z_marg_lamgam: $(Z_marg_lamgam)")
    println("use_repulsive: $(use_repulsive)")
    println("Z_thin: $(Z_thin)")
    flush(stdout)
  end

  @assert printFreq >= -1
  if printFreq == 0
    numPrints = 10
    printFreq = Int(ceil((nburn + nmcmc) / numPrints))
  end

  # Initialize tuners if needed
  if tuners == nothing
    t = Tuners(d.data.y, c.constants.K)
    TunersFS(t, init.theta, d.data.X)
  end

  function printMsg(iter::Int, msg::String)
    if printFreq > 0 && iter % printFreq == 0
      print(msg)
    end
  end

  # Loglike
  loglike = Float64[]

  # Instantiate (but not initialize) CPO stream
  if computeLPML
    cpoStream = MCMC.CPOstream{Float64}()
  end

  # DIC
  if computeDIC
    local tmp = DICparam(p=deepcopy(d.data.y),
                         mu=deepcopy(d.data.y),
                         sig=[zeros(Float64, d.data.N[i]) for i in 1:d.data.I],
                         y=deepcopy(d.data.y))
    dicStream = MCMC.DICstream{DICparam}(tmp)

    function updateParams(d::MCMC.DICstream{DICparam}, param::DICparam)
      d.paramSum.p += param.p
      d.paramSum.mu += param.mu
      d.paramSum.sig += param.sig
      d.paramSum.y += param.y

      return
    end

    function paramMeanCompute(d::MCMC.DICstream{DICparam})::DICparam
      return DICparam(d.paramSum.p / d.counter,
                      d.paramSum.mu / d.counter,
                      d.paramSum.sig / d.counter,
                      d.paramSum.y / d.counter)
    end

    function loglikeDIC(param::DICparam)::Float64
      ll = 0.0

      for i in 1:d.data.I
        for j in 1:d.data.J
          for n in 1:d.data.N[i]
            y_inj_is_missing = (d.data.m[i][n, j] == 1)

            # NOTE: Refer to `../compute_loglike.jl` for reasoning.
            if y_inj_is_missing

              # Compute p(m_inj | y_inj, theta) term.
              ll += log(param.p[i][n, j])
            end

            # Compute p(y_inj | theta) term.
            ll += logpdf(Normal(param.mu[i][n, j], param.sig[i][n]),
                         param.y[i][n, j])
          end
        end
      end

      return ll
    end

    function convertStateToDicParam(s::State)::DICparam
      p = [[prob_miss(s.y_imputed[i][n, j], c.constants.beta[:, i])
            for n in 1:d.data.N[i], j in 1:d.data.J] 
           for i in 1:d.data.I]

      mu = [[mus(i, n, j, s, c.constants, d.data)
             for n in 1:d.data.N[i], j in 1:d.data.J] 
            for i in 1:d.data.I]

      sig = [fill(sqrt(s.sig2[i]), d.data.N[i]) for i in 1:d.data.I]

      y = deepcopy(s.y_imputed)

      return DICparam(p, mu, sig, y)
    end
  end

  dden = Matrix{Vector{Float64}}[]

  function update!(s::StateFS, iter::Int, out)
    update_state_feature_select!(s, c, d, tuners, 
                                 ll=loglike, fix=fix,
                                 use_repulsive=use_repulsive,
                                 Z_marg_lamgam=Z_marg_lamgam,
                                 sb_ibp=sb_ibp, time_updates=time_updates,
                                 Z_thin=Z_thin)

    if computedden && iter > nburn && (iter - nburn) % thin_dden == 0
      append!(dden,
              [[datadensity(i, j, s.theta, c.constants, d.data)
                for i in 1:d.data.I, j in 1:d.data.J]])
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

      DICg = MCMC.DIC_gelman(sum(d.data.N) * loglike[(nburn+1):end])
      printMsg(iter, " -- DIC_gelman: $(DICg)")
    end

    printMsg(iter, "\n")
    flush(stdout)
  end

  if isinf(compute_loglike(init.theta, c.constants, d.data, normalize=true))
    println("Warning: Initial state yields likelihood of zero.")
    msg = """
    It is likely the case that the initialization of missing values
    is not consistent with the provided missing mechanism. The MCMC
    will almost certainly reject the initial values and sample new
    ones in its place.
    """
    println(join(split(msg), " "))
  else
    println("")
  end

  out, lastState = MCMC.gibbs(init, update!, monitors=monitors,
                              thins=thins, nmcmc=nmcmc, nburn=nburn,
                              printFreq=printFreq,
                              loglike=loglike, 
                              printlnAfterMsg=false)

  mega_out = Dict{Symbol, Any}(:samples => out,
                               :lastState => lastState,
                               :loglike => loglike,
                               :c => c,
                               :init => init)

  if computeDIC || computeLPML
    LPML = computeLPML ? MCMC.computeLPML(cpoStream) : NaN
    Dmean, pD = computeDIC ? MCMC.computeDIC(dicStream,
                                             loglikeDIC,
                                             paramMeanCompute,
                                             return_Dmean_pD=true) : (NaN, NaN)
    DICg = MCMC.DIC_gelman(sum(d.data.N) * loglike[(nburn+1):end])

    metrics = Dict(:LPML => LPML,
                   :DIC => Dmean + pD,
                   :Dmean => Dmean,
                   :pD => pD,
                   :DICg => DICg)
    println()
    println("metrics:")
    for (k, v) in metrics
      mega_out[k] = v
      println("$k => $v")
    end
  end
  flush(stdout)

  if computedden
    mega_out[:dden] = dden
  end

  mega_out[:nburn] = nburn
  mega_out[:Z_thin] = Z_thin
  mega_out[:m] = Matrix{Bool}.(d.data.m)

  return mega_out
end
