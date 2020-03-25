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
function outdir_suffix(pthin, batchprop, alpha, Nfac)
  return "pthin$(pthin)-batchprop$(batchprop)-alpha$(alpha)-N$(Nfac)" # TODO:N2000
end


# NOTE: modify
settings = let
  Z = Matrix{Bool}(readdlm(joinpath(@__DIR__, "Z.txt"), Int, comments=true))
  W = readdlm(joinpath(@__DIR__, "W.txt"), comments=true)
  Nfac = 2000  # for Nfac=2000, wasn't able to recover feature 1.
  N = [1, 1] *  Nfac
  [Dict(:simname => simname,
        :repfam_dist_scale => 0.1,  # doesn't work for kmeans init when phi=1.0
        :N => N,
        :Z => Z,
        :W => W,
        :thin_samps => 2,
        :nburn => 500,
        :nsamps => 1000,
        :Lmcmc => Dict(0 => 2, 1 => 2),
        :Kmcmc => 10,
        :pthin => pthin,
        :alpha => alpha,
        :batchprop => batchprop,
        :dataseed => 1,
        :mcmcseed => 1,
        :outdir_suffix => outdir_suffix(pthin, batchprop, alpha, Nfac))
   for pthin in [2]  # pthin{2,10} didn't work for Nfac=2000
   for batchprop in [.05]  # batchprop{.05} didn't work for Nfac=2000
   for alpha in [1.0]]  # alpha{10} didn't work for Nfac=2000
end
