# NOTE: Carefully review all the NOTES here.

include(joinpath(@__DIR__, "../../imports.jl"))
include(joinpath(@__DIR__, "../../simulatedata2.jl"))

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
  sim_z = CytofRepFAM.Model.sim_fn_exp_decay_generator(phi)
  use_repulsive = phi > 0

  function init_state_const_data(simdat; K, L)
    deltaz_prior = TruncatedNormal(1.0, 0.1, 0.0, Inf)
    d = CytofRepFAM.Model.Data(simdat[:y])
    c = CytofRepFAM.Model.defaultConstants(d, K, L,
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

    # NOTE: This is new. We want to spread weights evenly.
    c.eta_prior = Dict(z => Dirichlet(L[z], 1) for z in 0:1)

    # Initialize state
    # s = CytofRepFAM.Model.smartInit(c, d)  # mclust init
    s = CytofRepFAM.Model.smartInit(c, d, modelNames="kmeans",
                                    seed=settings[:mcmcseed],
                                    iterMax=30)  # kmeans init
    # s = CytofRepFAM.Model.genInitialState(c, d,  # NOTE: random inits
    #                                       sb_ibp=false,
    #                                       allow_repeated_Z_columns=false)

    t = CytofRepFAM.Model.Tuners(d.y, c.K)
    X = CytofRepFAM.Model.eye(Float64, d.I)

    cfs = CytofRepFAM.Model.ConstantsFS(c)

    # NOTE: These priors are important
    # similar to p ~ Beta(1, 9), weakly informative
    cfs.omega_prior = Normal(-3, 1.3)
    cfs.W_star_prior = Gamma(1.0, 2) # shape, scale

    dfs = CytofRepFAM.Model.DataFS(d, X)
    sfs = CytofRepFAM.Model.StateFS{Float64}(s, dfs)
    tfs = CytofRepFAM.Model.TunersFS(t, s, X)

    return Dict(:dfs => dfs, :cfs => cfs, :sfs => sfs, :tfs => tfs,
                :simdat => simdat)
  end

  @time simdat = simulatedata2(Z=settings[:Z],
                               N=settings[:N],
                               W=settings[:W],  # Categorical.
                               sig2=[.5, .5],
                               mus=Dict(0=>[-1.0], 1=>[1.0]),
                               skew=-0.9,
                               seed=settings[:dataseed],
                               propmissingscale=0, # NOTE: no missing data
                               sortLambda=true, 
                               eps_mus_dist=Uniform(-.3, .3));

  # Parameters to monitor
  monitor1 = [:theta__Z, :theta__v, :theta__alpha,
              :omega, :r, :theta__lam, :W_star, :theta__eta,
              :theta__W, :theta__delta, :theta__sig2]
  monitor2 = [:theta__y_imputed, :theta__gam]

  # MCMC Specs
  nsamps_to_thin(nsamps::Int, nmcmc::Int) = max(1, div(nmcmc, nsamps))
  nsamps = settings[:nsamps]  # Number of samples
  thin_samps = settings[:thin_samps]  # Factor to thin the primary parameters
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
  CytofRepFAM.Model.printConstants(config[:cfs])
  flush(stdout)

  # Fit model
  @time out = CytofRepFAM.Model.fit_fs_tp!(
    config[:sfs], config[:cfs], config[:dfs],
    nmcmc=mcmc_iter, nburn=nburn,
    batchprop=settings[:batchprop],
    prior_thin=settings[:pthin],
    temper=1.0, anneal=true, mb_update_burn_prop=0.7,
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
    thin_dden=nsamps_to_thin(200, mcmc_iter),
    time_updates=false,
    verbose=3)

  # Dump output
  BSON.bson("$(results_dir)/output.bson", out)
  BSON.bson("$(results_dir)/simdat.bson", Dict(:simdat => config[:simdat]))

  println("Completed!")
end  # simfn
