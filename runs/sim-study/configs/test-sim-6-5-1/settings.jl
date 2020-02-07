#=
using DelimiteFiles
=#

# NOTE: These simulation setting are identical to that of test-sims-6-5 EXCEPT
# that W_true is different. The W_true here has low abundance (< 10% each) for
# 2 pairs (i.e. 4) features. Thus, they should be harder to recover, and
# uncertainty in R_i should be greater.

# Simulation Name
simname = basename("$(@__DIR__)")
println(simname); flush(stdout)

# NOTE: write to scratchdir
function outdir_suffix(scale, kmcmc, seed_data)
  return "seed_$(seed_data)/scale_$(scale)/Kmcmc_$(kmcmc)"
end


# NOTE: modify
settings = let
  Z = Matrix{Bool}(readdlm(joinpath(@__DIR__, "Z.txt"), Int, comments=true))
  W = readdlm(joinpath(@__DIR__, "W.txt"), comments=true)
  [Dict(:simname => simname,
        :repfam_dist_scale => scale,
        :N => [1000, 1000],
        :Z => Z,
        :W => W,
        :thin_samps => 2,
        :nburn => 10000,
        :nsamps => 2000,
        # :nburn => 10, # NOTE: test
        # :nsamps => 20, # NOTE: test
        :Lmcmc => Dict(0 => 2, 1 => 2),
        :Kmcmc => kmcmc,
        :seed_data => seed_data,
        :seed_mcmc => 0,
        :outdir_suffix => outdir_suffix(scale, kmcmc, seed_data))
   for scale in [0, 1, 10]
   for kmcmc in [5, 6, 7, 8, 15]
   for seed_data in 1:3]
end
