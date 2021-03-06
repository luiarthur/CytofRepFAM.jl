#= Load these if running these tests in terminal
import Pkg; Pkg.activate("../")  # CytofRepFAM
using CytofRepFAM, Random, BSON, Test, Distributions
=#

function init_state_const_data(; N=[300, 200, 100], J=8, K=4,
                               L=Dict(0=>5, 1=>3), KMCMC=K*2,
                               sig2=fill(0.1, length(N)),
                               mus=Dict(0=>collect(range(-5, length=L[0], stop=-1)),
                                        1=>collect(range(1, length=L[1], stop=5))),
                               LMCMC=Dict(0=>5, 1=>5), seed=-1,
                               allow_repeated_Z_columns=true)

  if seed >= 0
    Random.seed!(seed)
  end

  I = length(N)
  simdat = CytofRepFAM.Model.genData(J, N, K, L, 
                                     mus=mus, sig2=sig2,
                                     sortLambda=true)
  d = CytofRepFAM.Model.Data(simdat[:y])
  simz = CytofRepFAM.Model.sim_fn_exp_decay_generator(1.0)
  deltaz_prior = TruncatedNormal(1.0, 0.1, 0.0, Inf)
  c = CytofRepFAM.Model.defaultConstants(d, KMCMC, LMCMC, 
                                         # sig2_prior=InverseGamma(11, 5),
                                         sig2_prior=InverseGamma(3, 2),
                                         delta0_prior=deltaz_prior,
                                         delta1_prior=deltaz_prior,
                                         alpha_prior=Gamma(0.1, 10.0),
                                         yQuantiles=[.01, .1, .25],
                                         pBounds=[.4, .8, .05],
                                         similarity_Z=simz)
  s = CytofRepFAM.Model.genInitialState(c, d,
                                        sb_ibp=false,
                                        allow_repeated_Z_columns=allow_repeated_Z_columns)
  t = CytofRepFAM.Model.Tuners(d.y, c.K)
  X = CytofRepFAM.Model.eye(Float64, d.I)
  println("mu*0 truth: $(simdat[:mus][0])")
  println("mu*1 truth: $(simdat[:mus][1])")
  println("sig2 truth: $(simdat[:sig2])")

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
                             computedden=true, seed=0)

  outdir = "results/repfam"
  mkpath(outdir)

  BSON.bson("$(outdir)/out_fs.bson", out)
  BSON.bson("$(outdir)/data_fs.bson", Dict(:simdat => config[:simdat]))
end

