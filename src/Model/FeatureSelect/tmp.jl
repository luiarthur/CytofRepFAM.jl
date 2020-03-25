function Z_changed(W, Z1, Z2)
  selected_features = findall(vec(sum(W, dims=1)) .> 0)
  return all(Z1[:, selected_features] .== Z2[:, selected_features])
end

function update_via_trained_prior!(sfs, dfs, cfs, tfs,
                                   batchprop::Float64, prior_thin::Int;
                                   fix, use_repulsive, Z_marg_lamgam, sb_ibp,
                                   time_updates, temper::Float64=1.0,
                                   verbose::Int=0)

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
  batchsizes = round.(Int, dfs.data.N * batchprop)
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
    update_state_feature_select!(sfs_mini, cfs, dfs_mini, tfs,
                                 ll=zeros(0),
                                 fix=vcat(fix,
                                          [:y_imputed,
                                          :v, :alpha, :omega,
                                          :r, :lam, :W_star, :gam,
                                          :eta, :delta, :sig2,
                                          :marg_loglike]),
                                 use_repulsive=use_repulsive,
                                 Z_marg_lamgam=Z_marg_lamgam,
                                 sb_ibp=sb_ibp, time_updates=time_updates,
                                 verbose=0)
  end

  # FIXME!
  # Compute metropolis ratio
  log_accept_ratio = if Z_changed(sfs.theta.W, sfs.theta.Z, sfs_mini.theta.Z)
    y_mega = [deepcopy(sfs.theta.y_imputed[i][idx_mega[i], :])
              for i in 1:dfs.data.I]
    c = cfs.constants
    ll_prop = compute_marg_loglike(sfs_mini.theta, c, d_mega, 1.0,
                                   y=y_mega, compute_prob_miss=false)
    ll_curr = compute_marg_loglike(sfs.theta, c, d_mega, 1.0,
                                   y=y_mega, compute_prob_miss=false)
    if verbose > 0
      print(" -- Z changed!")
      println("\n old Z:")
      Base.show(stdout, "text/plain", sfs.theta.Z)
      println("\n new Z:")
      Base.show(stdout, "text/plain", sfs_mini.theta.Z)
      print(" -- ll_prop: $(round(ll_prop, digits=2))")
      print(" -- ll_curr: $(round(ll_curr, digits=2))")
    end
    (ll_prop - ll_curr) / temper
  else
    if verbose > 0
      print(" -- Z unchanged!")
    end
    0.0
  end

  if verbose > 0
    print(" -- log_acceptance_ratio: $(round(log_accept_ratio, digits=2))")
  end
  if log_accept_ratio > log(rand())
    if verbose > 0
      print(" (accepted)")
    end
    sfs.theta.Z .= deepcopy(sfs.theta.Z)
    # curr_lam = deepcopy(sfs.theta.lam)
    # curr_gam = deepcopy(sfs.theta.gam)
    # curr_y = deepcopy(sfs.theta.y_imputed)

    # sfs.theta = sfs_mini.theta

    # sfs.theta.lam .= curr_lam
    # sfs.theta.gam .= curr_gam
    # sfs.theta.y_imputed .= curr_y

    # sfs.r .= deepcopy(sfs_mini.r)
    # sfs.W_star .= deepcopy(sfs_mini.W_star)
    # sfs.omega .= deepcopy(sfs_mini.omega)
  else
    if verbose > 0
      print(" (rejected)")
    end
  end

  update_state_feature_select!(sfs, cfs, dfs, tfs,
                               ll=zeros(0),
                               fix=vcat(fix, :Z),
                               use_repulsive=use_repulsive,
                               Z_marg_lamgam=Z_marg_lamgam,
                               sb_ibp=sb_ibp,
                               verbose=verbose,
                               time_updates=time_updates)


  # update_lam!(sfs.theta, cfs.constants, dfs.data)
  # update_gam!(sfs.theta, cfs.constants, dfs.data)
  # update_y_imputed!(sfs.theta, cfs.constants, dfs.data, tfs.tuners)

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
end