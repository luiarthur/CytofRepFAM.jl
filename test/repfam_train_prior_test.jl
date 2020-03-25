#= Load these if running these tests in terminal
import Pkg; Pkg.activate("..")  # CytofRepFAM
using CytofRepFAM, Random, BSON, Test, Distributions
include("repfamFS_tests.jl")
=#

printstyled("Test fitting repFAM via trained priors...\n", color=:yellow)
@testset "repFAM PT" begin
  config = init_state_const_data(N=[3,2]*500, L=Dict(0=>1, 1=>1),
                                 LMCMC=Dict(0=>3, 1=>3),
                                 mus=Dict(0 => [-2.0], 1 => [2.0]),
                                 sig2=fill(2.0, 2), seed=0)
  cfs = CytofRepFAM.Model.ConstantsFS(config[:c])
  dfs = CytofRepFAM.Model.DataFS(config[:d], config[:X])
  sfs = CytofRepFAM.Model.StateFS{Float64}(config[:s], dfs)
  tfs = CytofRepFAM.Model.TunersFS(config[:t], config[:s], config[:X])

  # CytofRepFAM.Model.printConstants(cfs)

  # Fit model.
  # For algo tests
  nmcmc = 200
  nburn = 100

  # For compile tests
  # nmcmc = 3
  # nburn = 3

  Nsum = sum(dfs.data.N)
  alpha = 1.0
  @time out = CytofRepFAM.Model.fit_fs_tp!(sfs, cfs, dfs,
                                           nmcmc=nmcmc, nburn=nburn,
                                           batchprop=0.05,
                                           prior_thin=4,
                                           # iMCMC settings
                                           temper=(alpha + Nsum) / alpha,
                                           anneal=false,
                                           mb_update_burn_prop=0.6,
                                           #
                                           Z_marg_lamgam=1.0,
                                           Z_marg_lamgam_decay_rate=100.0,
                                           # Z_marg_lamgam_min=0.05,
                                           Z_marg_lamgam_min=1.0,
                                           printFreq=1, seed=0,
                                           computedden=true,
                                           computeDIC=true,
                                           computeLPML=true,
                                           thin_dden=2,
                                           time_updates=false,
                                           verbose=3)
  println("Writing Output ...") 

  outdir = "results/repfam-tp"
  mkpath(outdir)

  BSON.bson("$(outdir)/out_fs_pt.bson", out)
  BSON.bson("$(outdir)/data_fs_pt.bson", Dict(:simdat => config[:simdat]))
end

#= Test
outdir = "results/repfam-tp"
out = BSON.load(joinpath(outdir, "out_fs_pt.bson"))
simdat = BSON.load(joinpath(outdir, "data_fs_pt.bson"))[:simdat]
samps = out[:samples][1]
extract(sym) = [s[sym] for s in samps]
Zs = extract(:theta__Z)
mean(Zs)
Ws = extract(:theta__W)
mean(Ws)
etas = extract(:theta__eta)

deltas = extract(:theta__delta)
acceptance_rate = length(unique(deltas)) / length(deltas)
mus0s = -[cumsum(d[0]) for d in deltas]
mus1s = [cumsum(d[1]) for d in deltas]
mus0 = hcat(mus0s...)'
mus1 = hcat(mus1s...)'
sig2s = extract(:theta__sig2)
sig2 = hcat(sig2s...)'

using RCall
@rimport rcommon
@rimport grDevices
@rimport graphics
begin
  nburn = out[:nburn]
  grDevices.pdf("$(outdir)/results.pdf");
  rcommon.plotPosts(mus0, cname="mus0");
  rcommon.plotPosts(mus1, cname="mus1");
  rcommon.plotPosts(sig2, cname="sig2");
  graphics.plot(out[:ll][end-nburn:end], main="loglike", type="l",
                xlab="iter", ylab="loglike");
  grDevices.dev_off();
end


include("../runs/PlotUtils/PlotUtils.jl")
PlotUtils.plot_dden(ddens=out[:dden],
                    etas=etas, Ws=Ws, Zs=Zs, sig2s=sig2s, deltas=deltas,
                    ygrid=out[:c].constants.y_grid,
                    imgdir=outdir, simdat=simdat)
=#
