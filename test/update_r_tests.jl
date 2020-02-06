#= Load these if running these tests in terminal
import Pkg; Pkg.activate("../")  # CytofRepFAM
using CytofRepFAM, Random, BSON, Test, Distributions
=#

@testset "time update of r" begin
  Random.seed!(0)

  config = init_state_const_data(N=[3, 2, 1] * 10000)

  cfs = CytofRepFAM.Model.ConstantsFS(config[:c])
  dfs = CytofRepFAM.Model.DataFS(config[:d], config[:X])
  sfs = CytofRepFAM.Model.StateFS{Float64}(config[:s], dfs)
  tfs = CytofRepFAM.Model.TunersFS(config[:t], config[:s], config[:X])

  println("N: $(dfs.data.N)")

  println("updating r ...")

  # Warmup
  CytofRepFAM.Model.update_r!(sfs, cfs, dfs)

  # time these:
  for i in 1:5
    @time CytofRepFAM.Model.update_r!(sfs, cfs, dfs)
  end
  println("checksum: $(sum(sfs.theta.lam[2]))")
end
