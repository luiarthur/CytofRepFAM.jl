function update_via_trained_prior!(sfs, dfs, cfs, tfs,
                                   batchprop::Float64,
                                   prior_thin::Int;
                                   batchsizes=nothing,
                                   fix, use_repulsive, Z_marg_lamgam, sb_ibp,
                                   time_updates, temper::Float64=1.0,
                                   verbose::Int=0,
                                   minibatch_update_all_params=false)

  # TODO:
  # Need to chane this to:
  # 1. Sample theta given minibatch
  #     - Note batchsize
  #     - Note thinning factor
  # 2. Accept theta (subset) with metropolis ratio 
  # 3. Sample lambda
  # 4. Sample gamma
  # 5. Sample missing data

  # Generate minibatch
  if batchsizes == nothing
    batchsizes = round.(Int, dfs.data.N * batchprop)
  end
  idx_mini, idx_mega, d_mini, d_mega = sample_minibatch(batchsizes, dfs.data.y)
  dfs_mini = DataFS(d_mini, dfs.X)
  sfs_mini = deepcopy(sfs)
  for i in 1:dfs.data.I
    sfs_mini.theta.lam[i] = deepcopy(sfs_mini.theta.lam[i][idx_mini[i]])
    sfs_mini.theta.gam[i] = deepcopy(sfs_mini.theta.gam[i][idx_mini[i], :])
    sfs_mini.theta.y_imputed[i] = deepcopy(
      sfs.theta.y_imputed[i][idx_mini[i], :])
  end

  # Sample theta given minibatch
  for _ in 1:prior_thin
    if minibatch_update_all_params
      update_state_feature_select!(sfs_mini, cfs, dfs_mini, tfs,
                                   ll=zeros(0),
                                   fix=vcat(fix, :y_imputed, :marg_loglike),
                                   use_repulsive=use_repulsive,
                                   Z_marg_lamgam=Z_marg_lamgam,
                                   sb_ibp=sb_ibp, time_updates=time_updates)
    else
      update_Z_marg_lamgam!(sfs_mini.theta, cfs.constants, dfs_mini.data,
                            sb_ibp, use_repulsive=use_repulsive)
    end
  end

  # FIXME!
  # Compute metropolis ratio
  log_accept_ratio = let
    y_mega = [deepcopy(sfs.theta.y_imputed[i][idx_mega[i], :])
              for i in 1:dfs.data.I]
    c = cfs.constants
    ll_prop = compute_marg_loglike(sfs_mini.theta, c, d_mega, 1.0,
                                   y=y_mega, compute_prob_miss=false)
    ll_curr = compute_marg_loglike(sfs.theta, c, d_mega, 1.0,
                                   y=y_mega, compute_prob_miss=false)
    if verbose > 0
      print(" -- ll_prop: $(round(ll_prop, digits=2))")
      print(" -- ll_curr: $(round(ll_curr, digits=2))")
    end
    (ll_prop - ll_curr) / temper
  end

  if verbose > 1
    println()
    println("current mu0: $(-cumsum(sfs.theta.delta[0]))")
    println("current mu1: $(cumsum(sfs.theta.delta[1]))")
    println("current sig2: $(sfs.theta.sig2)")
    println("current w: $(sfs.theta.W)")
    println("proposed mu0: $(-cumsum(sfs_mini.theta.delta[0]))")
    println("proposed mu1: $(cumsum(sfs_mini.theta.delta[1]))")
    println("proposed sig2: $(sfs_mini.theta.sig2)")
    println("proposed w: $(sfs_mini.theta.W)")
  end

  if verbose > 0
    print(" -- log_acceptance_ratio: $(round(log_accept_ratio, digits=2))")
  end

  # Whether or not to accept proposal
  accept_proposal = log_accept_ratio > log(rand())

  if accept_proposal
    if verbose > 0
      print(" (accepted)")

      # Print updated Z
      if verbose > 0.5
        ks_old = findall(vec(sum(sfs.r, dims=1) .> 0))
        ks_new = findall(vec(sum(sfs_mini.r, dims=1) .> 0))
        Z_old = sfs.theta.Z[:, ks_old]
        Z_new = sfs_mini.theta.Z[:, ks_new]
        if Z_old != Z_new
          println("\nOld Z:")
          Base.show(stdout, "text/plain", Z_old)
          println("\nNew Z:")
          Base.show(stdout, "text/plain", Z_new)
        end
      end
    end

    if minibatch_update_all_params
      curr_lam = deepcopy(sfs.theta.lam)
      curr_gam = deepcopy(sfs.theta.gam)
      curr_y = deepcopy(sfs.theta.y_imputed)

      sfs.theta = sfs_mini.theta

      sfs.theta.lam .= curr_lam
      sfs.theta.gam .= curr_gam
      sfs.theta.y_imputed .= curr_y

      sfs.r .= deepcopy(sfs_mini.r)
      sfs.W_star .= deepcopy(sfs_mini.W_star)
      sfs.omega .= deepcopy(sfs_mini.omega)
    else
      sfs.theta.Z .= deepcopy(sfs_mini.theta.Z)
    end
  else
    if verbose > 0
      print(" (rejected)")
    end
  end

  if minibatch_update_all_params
    update_lam!(sfs.theta, cfs.constants, dfs.data)
    update_gam!(sfs.theta, cfs.constants, dfs.data)
    update_y_imputed!(sfs.theta, cfs.constants, dfs.data, tfs.tuners)
  else
    update_state_feature_select!(sfs, cfs, dfs, tfs,
                                 ll=zeros(0),
                                 fix=vcat(fix, :Z),
                                 use_repulsive=use_repulsive,
                                 Z_marg_lamgam=Z_marg_lamgam,
                                 sb_ibp=sb_ibp, time_updates=time_updates)
  end

  if verbose > 1
    println("\nomega: $(sfs.omega)")
    println("\np: $(compute_p(dfs.X, sfs.omega))")
    println("W:")
    Base.show(stdout, "text/plain", sfs.theta.W)
    println()
    println("R: $(sum(sfs.r, dims=2))")
  end


  if size(unique(sfs.theta.Z, dims=2), 2) < cfs.constants.K && use_repulsive
    println("This Z has repeated columns after updating!")
    Base.show(stdout, "text/plain", sfs.theta.Z)
    println("\nFlipping $(cfs.constants.K) bits")
    sfs.theta.Z .= flip_bits(sfs.theta.Z, cfs.constants.K)
  end
end
