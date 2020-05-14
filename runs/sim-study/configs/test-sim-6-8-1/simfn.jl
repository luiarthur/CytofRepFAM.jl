# NOTE: Carefully review all the NOTES here.

include(joinpath(@__DIR__, "../../imports.jl"))
include(joinpath(@__DIR__, "simulatedata.jl"))
const rfam = CytofRepFAM.Model

function simfn(settings::Dict{Symbol, Any})
  println("pid: $(getpid())")
  println("Threads: $(Threads.nthreads())")
  println("pwd: $(pwd())")
  flush(stdout)

  println("settings:")
  for (k, v) in settings
    println("$k => $v")
  end

  # Results dir
  results_dir = settings[:results_dir]
  mkpath(results_dir)

  # Repulsive FAM penalty scale
  phi = settings[:repfam_dist_scale]
  sim_z = rfam.sim_fn_exp_decay_generator(phi)
  use_repulsive = phi > 0

  function init_state_const_data(simdat; K, L)
    deltaz_prior = TruncatedNormal(1.0, 0.1, 0.0, Inf)
    d = rfam.Data(simdat[:y])
    c = rfam.defaultConstants(d, K, L,
                              sig2_prior=InverseGamma(3, 1),
                              delta0_prior=deltaz_prior,
                              delta1_prior=deltaz_prior,
                              alpha_prior=Gamma(0.1, 10.0),
                              # NOTE:  These 3 items are important
                              yQuantiles=[.0, .25, .5], 
                              pBounds=[.05, .8, .05],
                              y_grid=collect(range(-10,
                                                   stop=4,
                                                   length=100)),
                              probFlip_Z=1.0,
                              similarity_Z=sim_z)

    # NOTE: To spread weights evenly.
    c.eta_prior = Dict(z => Dirichlet(L[z], 1) for z in 0:1)

    # Initialize state.
    # s = rfam.smartInit(c, d)  # mclust init
    states = [rfam.smartInit(c, d, modelNames="kmeans",
                             seed=settings[:mcmcseed] + s,
                             iterMax=100)  # kmeans init
              for s in 1:settings[:ntemps]]

    t = rfam.Tuners(d.y, c.K)
    X = rfam.eye(Float64, d.I)

    cfs = rfam.ConstantsFS(c)

    # NOTE: These priors are important
    # similar to p ~ Beta(1, 9), weakly informative
    cfs.omega_prior = Normal(-3, 1.3)
    cfs.W_star_prior = Gamma(1.0, 2) # shape, scale

    dfs = rfam.DataFS(d, X)
    sfss = [rfam.StateFS{Float64}(s, dfs) for s in states]
    tfs = rfam.TunersFS(t, states[1], X)

    return Dict(:dfs => dfs, :cfs => cfs, :sfss => sfss, :tfs => tfs,
                :simdat => simdat)
  end

  @time simdat = simulatedata(Z=settings[:Z],
                              N=settings[:N],
                              W=settings[:W],  # Categorical.
                              sig2=[.5, .5],
                              mus=Dict(0=>[-0.8], 1=>[1.3]),
                              skew=-0.9,
                              seed=settings[:dataseed],
                              propmissingscale=settings[:propmissingscale],
                              sortLambda=true, 
                              expressed_not_missing=true,
                              eps_mus_dist=Uniform(-.3, .3));

  # Parameters to monitor
  monitor1 = [:theta__Z, :theta__v, :theta__alpha,
              :omega, :r, :theta__lam, :W_star, :theta__eta,
              :theta__W, :theta__delta, :theta__sig2]
  monitor2 = [:theta__y_imputed, :theta__gam]

  # MCMC Specs
  nsamps_to_thin(nsamps::Int, nmcmc::Int) = max(1, div(nmcmc, nsamps))
  nsamps = settings[:nsamps]  # Number of samples
  thin_samps = settings[:thin_samps]  # Factor to thin the
                                      # primary parameters
  mcmc_iter = nsamps * thin_samps  # Number of MCMC iterations

  # LPML / DIC are computed based on `mcmc_iter` samples
  nburn = settings[:nburn]  # burn-in time

  # Configurations: priors, initial state, data, etc.
  config = init_state_const_data(simdat,
                                 K=settings[:Kmcmc],
                                 L=settings[:Lmcmc])

  # Print constants
  println("N: $(config[:dfs].data.N)")
  println("J: $(config[:dfs].data.J)")
  rfam.printConstants(config[:cfs])
  flush(stdout)

  # Fit model
  @time out = rfam.fit_fs_imcmc_pt!(
    config[:cfs], config[:dfs],
    inits=config[:sfss],
    nmcmc=mcmc_iter, nburn=nburn,
    batchprop=settings[:batchprop],
    prior_thin=settings[:pthin],
    tempers=settings[:temperatures],
    randpair=1.0, swap_freq=1.0,
    Z_marg_lamgam=1.0,
    Z_marg_lamgam_decay_rate=100.0,
    Z_marg_lamgam_min=1.0,
    printFreq=1,
    seed=settings[:mcmcseed],
    computedden=true,
    computeDIC=true,
    computeLPML=true,
    thins=[settings[:thin_samps],
           nsamps_to_thin(10, mcmc_iter)],
    ndden_samps=200,
    time_updates=false,
    verbose=3)

  # Dump output
  BSON.bson("$(results_dir)/output.bson", out)
  BSON.bson("$(results_dir)/simdat.bson",
            Dict(:simdat => config[:simdat]))

  println("Completed!")
end  # simfn
