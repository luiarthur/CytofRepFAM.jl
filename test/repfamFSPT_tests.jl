#= Load these if running these tests in terminal
import Pkg; Pkg.activate("..")  # CytofRepFAM
using CytofRepFAM, Random, BSON, Test, Distributions
using Distributed
include("repfamFS_tests.jl")
=#
using CytofRepFAM
using Distributed
using DelimitedFiles


printstyled("Test fitting repFAM on simulated data with PT...\n", color=:yellow)
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
  nmcmc = 1000
  nburn = 100
  maxcores = 20

  # For compile tests
  # nmcmc = 3
  # nburn = 3
  # maxcores = 2

  if length(workers()) != maxcores
    rmprocs(filter(w -> w > 1, workers()))
    addprocs(maxcores)
    @everywhere begin
      import Pkg; Pkg.activate("../")
      using CytofRepFAM
    end
  end

  # Temperatures
  # tempers = 1.025 .^ collect(0:(maxcores-1))
  # tempers = 1.02 .^ collect(0:(maxcores-1))
  # tempers = 1.1 .^ collect(0:(maxcores-1))
  # tempers = 1.01 .^ collect(0:(maxcores-1))  # swaps a lot
  # tempers = 2.0 .^ collect(0:(maxcores-1))
  # tempers = 1.1 .^ collect(0:(maxcores-1))
  # tempers = fill(1.0, maxcores)
  # tempers = 2.0 .^ ((collect(1:maxcores) .- 1) / (maxcores - 1))
  # tempers = 10.0 .^ ((collect(1:maxcores) .- 1) / (maxcores - 1))
  # tempers = 100.0 .^ ((collect(1:maxcores) .- 1) / (maxcores - 1))
  function gentempers(maxtemp, ntemps; degree=1)
    i = collect(1:ntemps)
    powers = (i.^ degree .- 1) / (ntemps ^ degree - 1)
    return maxtemp .^ powers
  end
  tempers = gentempers(50, maxcores, degree=2)

  @time out = CytofRepFAM.Model.fit_fs_pt!(cfs, dfs,
                                           tempers=tempers,
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

  @time out = CytofRepFAM.Model.fit_fs_pt!(cfs, dfs,
                                           # inits=[deepcopy(sfs) for _ in tempers],
                                           swap_freq=.5,
                                           tempers=tempers,
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
  # rmprocs(filter(w -> w > 1, workers()))

  println("Writing Output ...") 

  outdir = "results/repfam-pt"
  mkpath(outdir)

  open("$(outdir)/tempers.txt", "w") do io
    writedlm(io, tempers)
  end

  for t in 1:length(out[:lls])
    open(joinpath(outdir, "ll_$(t).txt"), "w") do io
      writedlm(io, out[:lls][t], ',')
    end
  end

  open(joinpath(outdir, "swapcounts.txt"), "w") do io
    writedlm(io, out[:swapcounts], ',')
  end

  open(joinpath(outdir, "paircounts.txt"), "w") do io
    writedlm(io, out[:paircounts], ',')
  end

  BSON.bson("$(outdir)/out_fs_pt.bson", out)
  BSON.bson("$(outdir)/data_fs_pt.bson", Dict(:simdat => config[:simdat]))
end
