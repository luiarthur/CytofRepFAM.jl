#= Load these if running these tests in terminal
import Pkg; Pkg.activate("../")  # CytofRepFAM
using CytofRepFAM, Random, BSON, Test, Distributions
=#

printstyled("Test fitting repFAM on simulated data with PT...\n", color=:yellow)
@testset "repFAM PT" begin
  config = init_state_const_data() 
  cfs = CytofRepFAM.Model.ConstantsFS(config[:c])
  dfs = CytofRepFAM.Model.DataFS(config[:d], config[:X])
  sfs = CytofRepFAM.Model.StateFS{Float64}(config[:s], dfs)
  tfs = CytofRepFAM.Model.TunersFS(config[:t], config[:s], config[:X])

  # CytofRepFAM.Model.printConstants(cfs)

  # Fit model.
  # For algo tests
  # nmcmc = 100
  # nburn = 1200

  # For compile tests
  nmcmc = 3
  nburn = 5
  out = CytofRepFAM.Model.fit_fs_pt!(sfs, cfs, dfs, tuners=tfs,
                                     tempers=[1.0, 2.0], ncores=2,
                                     nmcmc=nmcmc, nburn=nburn,
                                     printFreq=1, time_updates=true,
                                     computeDIC=true, computeLPML=true,
                                     Z_thin=5, seed=0)

  # outdir = "results/repfam"
  # mkpath(outdir)

  # BSON.bson("$(outdir)/out_fs_pt.bson", out)
  # BSON.bson("$(outdir)/data_fs_pt.bson", Dict(:simdat => config[:simdat]))
end

