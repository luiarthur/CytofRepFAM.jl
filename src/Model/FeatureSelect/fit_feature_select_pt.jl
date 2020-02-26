nsamps_to_thin(nsamps::Int, nmcmc::Int) = max(1, div(nmcmc, nsamps))

monitor1 = [:theta__Z, :theta__v, :theta__alpha,
            :omega, :r, :theta__lam, :W_star, :theta__eta,
            :theta__W, :theta__delta, :theta__sig2]

monitor2 = [:theta__y_imputed, :theta__gam]

# function fit_fs_pt!(init::StateFS, c::ConstantsFS, d::DataFS;
#                     tempers::Vector{Float64}, ncores::Int,
#                     nmcmc::Int, nburn::Int, 
#                     swap_freq::Int=10,
#                     remove_current_workers::Bool=true,
#                     tuners::Union{Nothing, TunersFS}=nothing,
#                     monitors=[monitor1, monitor2],
#                     fix::Vector{Symbol}=Symbol[],
#                     thins::Vector{Int}=[2, nsamps_to_thin(10, nmcmc)],
#                     thin_dden::Int=1,
#                     printFreq::Int=0, 
#                     computeDIC::Bool=false, computeLPML::Bool=false,
#                     computedden::Bool=false,
#                     sb_ibp::Bool=false,
#                     use_repulsive::Bool=true, joint_update_Z::Bool=true,
#                     verbose::Int=1, time_updates::Bool=false, Z_thin::Int=0,
#                     seed::Int=-1)
# 
#   # Remove current workers
#   if remove_current_workers
#     rmprocs(filter(w -> w > 1, workers()))
#   end
# 
#   # Spawn `ncores` new processes
#   addprocs(ncores)
# 
#   # TODO: Begin
#   @everywhere function update!(s, c, d, t)
#     loglike = Float64[]
#     update_state_feature_select!(s, c, d, t, loglike,
#                                  Symbol[], true, false, false,
#                                  time_updates=false)
#   end
# 
#   println("HERE")
#   _ = pmap(update!,
#            [deepcopy(init) for _ in tempers],
#            [deepcopy(c) for _ in tempers],
#            [deepcopy(d) for _ in tempers],
#            [deepcopy(tuners) for _ in tempers])
#   println("THERE")
#   # TODO: End
# 
# 
#   # Remove current workers
#   if remove_current_workers
#     rmprocs(filter(w -> w > 1, workers()))
#   end
# end
