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
function outdir_suffix(dataseed, mcmcseed, phi, Zind)
  return "dataseed$(dataseed)-mcmcseed$(mcmcseed)-phi$(phi)-Zind$(Zind)"
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
  [Dict(:simname => simname,
        :repfam_dist_scale => phi,
        :N => N,
        :Z => Z[Zind],
        :W => W,
        :thin_samps => 1,
        :nburn => 6000,
        :nsamps => 3000,
        :Lmcmc => Dict(0 => 3, 1 => 3),
        :Kmcmc => 10,
        :pthin => pthin,
        :batchprop => batchprop,
        :dataseed => dataseed,
        :mcmcseed => mcmcseed,
        :outdir_suffix => outdir_suffix(dataseed, mcmcseed, phi, Zind))
   for Zind in 1:3
   for phi in [0.0, 1.0, 10.0]
   for dataseed in 1:1
   for mcmcseed in 1:3]
end
