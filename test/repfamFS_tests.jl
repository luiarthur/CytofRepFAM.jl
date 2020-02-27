#= Load these if running these tests in terminal
import Pkg; Pkg.activate("../")  # CytofRepFAM
using CytofRepFAM, Random, BSON, Test, Distributions
=#

function init_state_const_data(; N=[300, 200, 100], J=8, K=4,
                               L=Dict(0=>5, 1=>3), seed=-1)

  if seed >= 0
    Random.seed!(seed)
  end

  I = length(N)
  simdat = CytofRepFAM.Model.genData(J, N, K, L, sortLambda=true)
  d = CytofRepFAM.Model.Data(simdat[:y])
  c = CytofRepFAM.Model.defaultConstants(d, K * 2, Dict(0=>5, 1=>5))
  s = CytofRepFAM.Model.genInitialState(c, d)
  t = CytofRepFAM.Model.Tuners(d.y, c.K)
  X = CytofRepFAM.Model.eye(Float64, d.I)

  _ = CytofRepFAM.Model.compute_marg_loglike(s, c, d, 1.0)
  @time ll = CytofRepFAM.Model.compute_marg_loglike(s, c, d, 1.0)
  println("Loglike with temper = 1.0, marginalized over lam, gam: $(ll)")

  @time ll = CytofRepFAM.Model.compute_marg_loglike(s, c, d, 2.0)
  println("Loglike with temper = 2.0, marginalized over lam, gam: $(ll)")


  return Dict(:d => d, :c => c, :s => s, :t => t, :X => X,
              :simdat => simdat)
end


printstyled("Test fitting repFAM on simulated data ...\n", color=:yellow)
@testset "update W*, r, omega" begin
  config = init_state_const_data(seed=0)
  cfs = CytofRepFAM.Model.ConstantsFS(config[:c])
  dfs = CytofRepFAM.Model.DataFS(config[:d], config[:X])
  sfs = CytofRepFAM.Model.StateFS{Float64}(config[:s], dfs)
  tfs = CytofRepFAM.Model.TunersFS(config[:t], config[:s], config[:X])

  # Do one update for W_star, r, omega.
  println("r init: $(sfs.r)")
  println("W_star init: $(sfs.W_star)")
  println("omega init: $(sfs.omega)")

  CytofRepFAM.Model.update_W_star!(sfs, cfs, dfs, tfs)
  CytofRepFAM.Model.update_r!(sfs, cfs, dfs)
  CytofRepFAM.Model.update_omega!(sfs, cfs, dfs, tfs)

  println("W_star after: $(sfs.W_star)")
  println("W after: $(sfs.theta.W)")
  println("r after: $(sfs.r)")
  println("omega after: $(sfs.omega)")

  CytofRepFAM.Model.printConstants(cfs)

  # Fit model.
  # For algo tests
  # nmcmc = 100
  # nburn = 1200
  # For compile tests
  nmcmc = 3
  nburn = 5
  out = CytofRepFAM.Model.fit_fs!(sfs, cfs, dfs, tuners=tfs,
                             nmcmc=nmcmc, nburn=nburn,
                             printFreq=1, time_updates=true,
                             computeDIC=true, computeLPML=true, Z_thin=1,
                             seed=0)

  outdir = "results/repfam"
  mkpath(outdir)

  BSON.bson("$(outdir)/out_fs.bson", out)
  BSON.bson("$(outdir)/data_fs.bson", Dict(:simdat => config[:simdat]))
end

