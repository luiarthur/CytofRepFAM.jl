function fit_fs_pt!(init::StateFS, cfs::ConstantsFS, dfs::DataFS, tfs::TunersFS;
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

  # Create distributed arrays for:
  # - loglike, ConstantsFS, StatesFS, TunersFS
  lls = distribute([Float64[0.0] for _ in tempers])
  cs = distribute([let 
                      ct = deepcopy(cfs)
                      ct.constants.temper = tau
                      ct
                    end for tau in tempers])
  states = distribute([deepcopy(init) for _ in tempers])
  tuners = distribute([deepcopy(tfs) for _ in tempers])

  println("About to start parallel chains ...")
  for iter in 1:(nburn+nmcmc)
    # Perform updates in parallel
    @sync @distributed for t in 1:num_tempers
      update_state_feature_select!(states[t], cs[t], dfs, tuners[t],
                                   ll=lls[t], fix=fix,
                                   use_repulsive=use_repulsive,
                                   Z_marg_lamgam=Z_marg_lamgam,
                                   sb_ibp=sb_ibp, time_updates=time_updates[t])
      println(lls[t])
    end

    if iter % swap_freq == 0
      println("TODO: PERFORM SWAP.")
    end
  end
  println("After running parallel chains ...")
end
