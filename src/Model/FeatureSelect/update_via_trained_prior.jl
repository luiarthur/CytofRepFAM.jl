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
    sfs_mini.theta.lam = [sfs.theta.lam[i][idx_mini[i]] for i in 1:dfs.data.I]
    sfs_mini.theta.gam = [sfs.theta.gam[i][idx_mini[i], :] for i in 1:dfs.data.I]

    # Sample theta given minibatch
    for _ in 1:prior_thin
      update_state_feature_select!(sfs_mini, cfs, dfs_mini, tfs,
                                   ll=zeros(), fix=vcat(fix, :y_imputed),
                                   use_repulsive=use_repulsive,
                                   Z_marg_lamgam=Z_marg_lamgam,
                                   sb_ibp=sb_ibp, time_updates=time_updates)
    end

    # Compute metropolis ratio
    log_accept_ratio = let
      dfs_mega = DataFS(d_mega, dfs.X)
      ll_prop = compute_marg_loglike(sfs_mini.theta, cfs, dfs_mega.data, 1.0)
      ll_curr = compute_marg_loglike(sfs.theta, cfs, dfs_mega.data, 1.0)
      ll_prop - ll_curr
    end

    if log_accept_ratio > log(rand())
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
    end

    update_lam!(sfs.theta, cfs.constants, dfs.data)
    update_gam!(sfs.theta, cfs.constants, dfs.data)
    update_y_imputed!(sfs.theta, cfs.constants, dfs.data, dft.tuners)
end
