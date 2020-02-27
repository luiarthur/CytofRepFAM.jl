nsamps_to_thin(nsamps::Int, nmcmc::Int) = max(1, div(nmcmc, nsamps))

monitor1 = [:theta__Z, :theta__v, :theta__alpha,
            :omega, :r, :theta__lam, :W_star, :theta__eta,
            :theta__W, :theta__delta, :theta__sig2]

monitor2 = [:theta__y_imputed, :theta__gam]

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
                    verbose::Int=1, time_updates::Bool=false, Z_thin::Int=1,
                    seed::Int=-1)

  println("HERE")

  # Number of temperatures
  num_tempers = length(tempers)

  # Create arrays
  lls = distribute([Float64[0.0] for _ in tempers])
  cs = distribute([let 
                      ct = deepcopy(cfs)
                      ct.constants.temper = tau
                      ct
                    end for tau in tempers])
  states = distribute([deepcopy(init) for _ in tempers])
  tuners = distribute([deepcopy(tfs) for _ in tempers])

  # Perform updates in parallel
  # pmap(states, cs, tunerss, lls) do s, c, tuner, ll
  #   println(time())
  #   update_state_feature_select!(s, c, d, tuner,
  #                                ll=ll, fix=fix,
  #                                use_repulsive=use_repulsive,
  #                                Z_marg_lamgam=Z_marg_lamgam,
  #                                sb_ibp=sb_ibp, time_updates=false)
  # end

  println("HERE2")
  @sync @distributed for t in 1:num_tempers
    println("Hey")
    update_state_feature_select!(states[t], cs[t], dfs, tuners[t],
                                 ll=lls[t], fix=Symbol[],
                                 use_repulsive=true,
                                 Z_marg_lamgam=true,
                                 sb_ibp=false, time_updates=false)
  end

  println("THERE")
end
