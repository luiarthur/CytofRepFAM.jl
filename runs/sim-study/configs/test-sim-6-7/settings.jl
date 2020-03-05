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
function outdir_suffix(maxtemp, ntemps, degree)
  return "maxtemp$(maxtemp)-ntempts$(ntemps)-degree$(degree)"
end


# NOTE: modify
settings = let
  Z = Matrix{Bool}(readdlm(joinpath(@__DIR__, "Z.txt"), Int, comments=true))
  W = readdlm(joinpath(@__DIR__, "W.txt"), comments=true)
  [Dict(:simname => simname,
        :repfam_dist_scale => 1.0,
        :N => [2000, 2000],
        :Z => Z,
        :W => W,
        :thin_samps => 2,
        :nburn => 10000,
        :nsamps => 2000,
        :Lmcmc => Dict(0 => 2, 1 => 2),
        :Kmcmc => 10,
        :maxtemp => maxtemp,
        :ntemps => ntemps,
        :ncores => ntemps,
        :degree => degree,
        :dataseed => 1,
        :mcmcseed => 1,
        :outdir_suffix => outdir_suffix(maxtemp, ntemps, degree))
   for maxtemp in [1000]
   for ntemps in [20]
   for degree in [4]]
end
