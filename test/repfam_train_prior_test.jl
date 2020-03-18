#= Load these if running these tests in terminal
import Pkg; Pkg.activate("..")  # CytofRepFAM
using CytofRepFAM, Random, BSON, Test, Distributions
include("repfamFS_tests.jl")
=#

printstyled("Test fitting repFAM via trained priors...\n", color=:yellow)
@testset "repFAM PT" begin
  config = init_state_const_data(N=[3,2]*100, L=Dict(0=>1, 1=>1),
                                 LMCMC=Dict(0=>2, 1=>2),
                                 mus=Dict(0 => [-2.0], 1 => [2.0]),
                                 sig2=fill(0.7, 2), seed=0)
  cfs = CytofRepFAM.Model.ConstantsFS(config[:c])
  dfs = CytofRepFAM.Model.DataFS(config[:d], config[:X])
  sfs = CytofRepFAM.Model.StateFS{Float64}(config[:s], dfs)
  tfs = CytofRepFAM.Model.TunersFS(config[:t], config[:s], config[:X])

  # CytofRepFAM.Model.printConstants(cfs)

  # Fit model.
  # For algo tests
  # nmcmc = 1000
  # nburn = 100
  # maxcores = 20

  # For compile tests
  nmcmc = 3
  nburn = 3

  @time out = CytofRepFAM.Model.fit_fs_tp!(sfs, cfs, dfs,
                                           nmcmc=3, nburn=3,
                                           Z_marg_lamgam=0.5,
                                           Z_marg_lamgam_decay_rate=10.0,
                                           Z_marg_lamgam_min=0.05,
                                           randpair=0.8,
                                           printFreq=1, seed=0,
                                           computedden=true,
                                           computeDIC=true,
                                           computeLPML=true,
                                           time_updates=false,
                                           verbose=1)

  @time out = CytofRepFAM.Model.fit_fs_tp!(cfs, dfs,
                                           nmcmc=nmcmc, nburn=nburn,
                                           Z_marg_lamgam=0.5,
                                           Z_marg_lamgam_decay_rate=100.0,
                                           Z_marg_lamgam_min=0.05,
                                           randpair=0.33,
                                           printFreq=1, seed=0,
                                           computedden=true,
                                           computeDIC=true,
                                           computeLPML=true,
                                           time_updates=false,
                                           verbose=3)
  println("Writing Output ...") 

  outdir = "results/repfam-tp"
  mkpath(outdir)

  BSON.bson("$(outdir)/out_fs_pt.bson", out)
  BSON.bson("$(outdir)/data_fs_pt.bson", Dict(:simdat => config[:simdat]))
end
