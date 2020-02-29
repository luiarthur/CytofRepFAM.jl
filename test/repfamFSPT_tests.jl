#= Load these if running these tests in terminal
import Pkg; Pkg.activate("..")  # CytofRepFAM
using CytofRepFAM, Random, BSON, Test, Distributions
using Distributed
include("repfamFS_tests.jl")
=#
using Distributed
using DelimitedFiles

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
  nmcmc = 200
  nburn = 200

  maxcores = 8
  rmprocs(filter(w -> w > 1, workers()))
  addprocs(maxcores)
  @everywhere begin
    import Pkg; Pkg.activate("../")
    using CytofRepFAM
  end

  # FIXME: fit_fs_pt.jl
  # tempers = 1.025 .^ collect(0:(maxcores-1))
  # tempers = 1.02 .^ collect(0:(maxcores-1))
  # tempers = 1.1 .^ collect(0:(maxcores-1))
  # tempers = 1.01 .^ collect(0:(maxcores-1))  # swaps a lot
  # tempers = 2.0 .^ collect(0:(maxcores-1))
  tempers = 1.1 .^ collect(0:(maxcores-1))
  # tempers = 1.1 .^ collect(0:(maxcores-1))
  # tempers = fill(1.0, maxcores)
  out = CytofRepFAM.Model.fit_fs_pt!(sfs, cfs, dfs, tfs,
                                     tempers=tempers,
                                     ncores=maxcores,
                                     swap_freq=1,
                                     nmcmc=2, nburn=2,
                                     printFreq=1, seed=0)
  @time out = CytofRepFAM.Model.fit_fs_pt!(sfs, cfs, dfs, tfs,
                                     tempers=tempers,
                                     swap_freq=2,
                                     ncores=maxcores,
                                     nmcmc=nmcmc, nburn=nburn,
                                     Z_marg_lamgam=false,
                                     printFreq=1, seed=0, verbose=2)
  rmprocs(filter(w -> w > 1, workers()))

  println("Writing Output ...") 

  outdir = "results/repfam-pt"
  mkpath(outdir)

  for t in 1:length(out[:lls])
    open(joinpath(outdir, "ll_$(t).txt"), "w") do io
      writedlm(io, out[:lls][t], ',')
    end
  end

  BSON.bson("$(outdir)/out_fs_pt.bson", out)
  BSON.bson("$(outdir)/data_fs_pt.bson", Dict(:simdat => config[:simdat]))
end
