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
function outdir_suffix(dataseed, mcmcseed, phi)
  return "dataseed$(dataseed)-mcmcseed$(mcmcseed)-phi$(phi)"
end

# NOTE: modify
settings = let
  Z = Matrix{Bool}(readdlm(joinpath(@__DIR__, "Z.txt"), Int, comments=true))
  W = readdlm(joinpath(@__DIR__, "W.txt"), comments=true)
  Nfac = 2000  # for Nfac=2000, wasn't able to recover feature 1.
  N = [1, 1] *  Nfac
  # NOTE: these are actually quite important.
  # Need pthin to be at least 4, but batchprop can be small, if N is
  # sufficiently large.
  pthin = 5
  batchprop = 0.05
  [Dict(:simname => simname,
        :repfam_dist_scale => phi,  # doesn't work for kmeans init when phi=1.0
        :N => N,
        :Z => Z,
        :W => W,
        :thin_samps => 1,
        :nburn => 3000,
        :nsamps => 3000,
        :Lmcmc => Dict(0 => 2, 1 => 2),
        :Kmcmc => 10,
        :pthin => pthin,
        :batchprop => batchprop,
        :dataseed => dataseed,
        :mcmcseed => mcmcseed,
        :outdir_suffix => outdir_suffix(dataseed, mcmcseed, phi))
   for phi in [0.0, 1.0, 10.]  # pthin{2,10} didn't work for Nfac=2000
   for dataseed in 1:3  # batchprop{.05} didn't work for Nfac=2000
   for mcmcseed in 1:3]  # alpha{10} didn't work for Nfac=2000
end
