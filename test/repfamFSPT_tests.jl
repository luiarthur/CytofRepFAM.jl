#= Load these if running these tests in terminal
import Pkg; Pkg.activate("..")  # CytofRepFAM
using CytofRepFAM, Random, BSON, Test, Distributions
using Distributed
include("repfamFS_tests.jl")
=#
using Distributed

printstyled("Test fitting repFAM on simulated data with PT...\n", color=:yellow)
@testset "repFAM PT" begin
  config = init_state_const_data(N=[3,1,2]*100, seed=0)
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
  nmcmc = 50
  nburn = 50

  maxcores = 4
  rmprocs(filter(w -> w > 1, workers()))
  addprocs(maxcores)
  @everywhere begin
    import Pkg; Pkg.activate("../")
    using CytofRepFAM
  end

  # FIXME: fit_fs_pt.jl
  tempers = 2.0 .^ collect(0:(maxcores-1))
  out = CytofRepFAM.Model.fit_fs_pt!(sfs, cfs, dfs, tfs,
                                     tempers=tempers,
                                     ncores=maxcores,
                                     nmcmc=2, nburn=2,
                                     printFreq=1, seed=0)
  @time out = CytofRepFAM.Model.fit_fs_pt!(sfs, cfs, dfs, tfs,
                                     tempers=tempers,
                                     swap_freq=5,
                                     ncores=maxcores,
                                     nmcmc=nmcmc, nburn=nburn,
                                     printFreq=1, seed=0, verbose=2)
  rmprocs(filter(w -> w > 1, workers()))

  # outdir = "results/repfam"
  # mkpath(outdir)

  # BSON.bson("$(outdir)/out_fs_pt.bson", out)
  # BSON.bson("$(outdir)/data_fs_pt.bson", Dict(:simdat => config[:simdat]))
end
