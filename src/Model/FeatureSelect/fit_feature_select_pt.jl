function swap!(i, j, x)
  tmp = x[i]
  x[i] = x[j]
  x[j] = tmp
  return x
end

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
  lls = [Float64[0.0] for _ in tempers]
  cs = [let 
          ct = deepcopy(cfs)
          ct.constants.temper = tau
          ct
        end for tau in tempers]
  states = [deepcopy(init) for _ in tempers]
  tuners = [deepcopy(tfs) for _ in tempers]

  println("About to start parallel chains ...")
  for iter in 1:(nburn+nmcmc)
    # Perform updates in parallel
    @time out = pmap(states, cs, tuners, lls, time_updates) do s, c, t, ll, tu
      @time for _ in 1:swap_freq
        update_state_feature_select!(s, c, dfs, t,
                                     ll=ll, fix=fix,
                                     use_repulsive=use_repulsive,
                                     # TODO: make this more explicit, instead
                                     # of random. See `update_Z_v2`.
                                     Z_marg_lamgam=Z_marg_lamgam,
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

    if iter % swap_freq == 0
      println("TODO: PERFORM SWAP.")
    end
  end
  println("After running parallel chains ...")
end
