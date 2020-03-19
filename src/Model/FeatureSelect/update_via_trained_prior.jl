function update_via_trained_prior!(sfs, dfs, cfs, tfs,
                                   batchprop::Float64, prior_thin::Int;
                                   fix, use_repulsive, Z_marg_lamgam, sb_ibp,
                                   time_updates)

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
                                   ll=zeros(0), fix=vcat(fix, :y_imputed),
                                   use_repulsive=use_repulsive,
                                   Z_marg_lamgam=Z_marg_lamgam,
                                   sb_ibp=sb_ibp, time_updates=time_updates)
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
      print(" -- ll_prop: $(round(ll_prop, digits=2))")
      print(" -- ll_curr: $(round(ll_curr, digits=2))")
      ll_prop - ll_curr
    end

    print("-- log_acceptance_ratio: $(round(log_accept_ratio, digits=2))")
    if log_accept_ratio > log(rand())
      print(" (accepted)")
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
      print(" (rejected)")
    end

    update_lam!(sfs.theta, cfs.constants, dfs.data)
    update_gam!(sfs.theta, cfs.constants, dfs.data)
    update_y_imputed!(sfs.theta, cfs.constants, dfs.data, tfs.tuners)
end
