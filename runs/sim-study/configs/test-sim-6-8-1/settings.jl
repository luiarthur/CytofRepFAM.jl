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
function outdir_suffix(pmiss, phi, zind)
  return "pmiss$(pmiss)-phi$(phi)-zind$(zind)"
end

# NOTE: modify
settings = let
  Z1 = Matrix{Bool}(readdlm(joinpath(@__DIR__, "Z1.txt"), Int, comments=true))
  Z2 = Matrix{Bool}(readdlm(joinpath(@__DIR__, "Z2.txt"), Int, comments=true))
  Z3 = Matrix{Bool}(readdlm(joinpath(@__DIR__, "Z3.txt"), Int, comments=true))
  Z = [Z1, Z2, Z3]
  W = readdlm(joinpath(@__DIR__, "W.txt"), comments=true)
  Nfac = 2000
  N = [1, 1] *  Nfac
  pthin = 5
  batchprop = 0.05
  dataseed = 1
  mcmcseed = 1
  temperatures = [1.0, 1.003, 1.006, 1.01]

  [Dict(:simname => simname,
        :repfam_dist_scale => phi,
        :N => N,
        :Z => Z[zind],
        :W => W,
        :thin_samps => 1,
        :nburn => 6000,
        :nsamps => 3000,
        :Lmcmc => Dict(0 => 3, 1 => 3),
        :Kmcmc => 10,
        :pthin => pthin,
        :temperatures => temperatures,
        :ntemps => length(temperatures),
        :propmissingscale => propmissingscale,
        :batchprop => batchprop,
        :dataseed => dataseed,
        :mcmcseed => mcmcseed,
        :outdir_suffix => outdir_suffix(propmissingscale, phi, zind))
   for propmissingscale in [0.6, 0.0]
   for phi in [1.0, 0.0, 10.0, 25.0]
   for zind in 1:3]
end
