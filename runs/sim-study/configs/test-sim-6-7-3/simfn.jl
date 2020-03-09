include(joinpath(@__DIR__, "../../imports.jl"))
include(joinpath(@__DIR__, "../../simulatedata.jl"))

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
                                      sig2_prior=InverseGamma(11, 5),
                                      delta0_prior=deltaz_prior,
                                      delta1_prior=deltaz_prior,
                                      alpha_prior=Gamma(0.1, 10.0),
                                      yQuantiles=[.0, .25, .5], 
                                      pBounds=[.05, .8, .05],
                                      probFlip_Z=1.0,
                                      similarity_Z=sim_z)

    # NOTE: This is new. We want to spread weights evenly.
    c.eta_prior = Dict(z => Dirichlet(L[z], 5) for z in 0:1)

    # s = CytofRepFAM.Model.genInitialState(c, d)
    s = CytofRepFAM.Model.smartInit(c, d)
    t = CytofRepFAM.Model.Tuners(d.y, c.K)
    X = CytofRepFAM.Model.eye(Float64, d.I)

    cfs = CytofRepFAM.Model.ConstantsFS(c)
    dfs = CytofRepFAM.Model.DataFS(d, X)
    sfs = CytofRepFAM.Model.StateFS{Float64}(s, dfs)
    tfs = CytofRepFAM.Model.TunersFS(t, s, X)

    return Dict(:dfs => dfs, :cfs => cfs, :sfs => sfs, :tfs => tfs,
                :simdat => simdat)
  end

  @time simdat = simulatedata1(Z=settings[:Z],
                               N=settings[:N],
                               W=settings[:W],  # Categorical.
                               sig2=[.5, .5],
                               seed=settings[:dataseed],
                               propmissingscale=.6,
                               sortLambda=true, 
                               eps_mus_dist=Uniform(-.5, .5));

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

  # Temperatures
  tempers = let
    nrep = 2
    println("nrep: $(nrep)")
    _tempers = CytofRepFAM.MCMC.WSPT.gentempers(settings[:maxtemp],
                                                settings[:ntemps] - nrep,
                                                degree=settings[:degree])
    append!(_tempers, ones(nrep))
    sort(_tempers)
  end

  println("tempers:")
  println(tempers)

  # Fit model
  @time out = CytofRepFAM.Model.fit_fs_pt!(
    config[:cfs], config[:dfs],
    tempers=tempers,
    nmcmc=mcmc_iter, nburn=nburn,
    inits=[deepcopy(config[:sfs]) for _ in tempers],
    swap_freq=.5,
    save_all_states=true,
    Z_marg_lamgam=0.1,
    Z_marg_lamgam_decay_rate=100.0,
    Z_marg_lamgam_min=0.05,
    randpair=0.3,
    thins=[thin_samps, nsamps_to_thin(10, mcmc_iter)],
    monitors=[monitor1, monitor2],
    computedden=true,
    thin_dden=nsamps_to_thin(200, mcmc_iter),
    printFreq=10, time_updates=false,
    computeDIC=true, computeLPML=true,
    use_repulsive=use_repulsive,
    verbose=3,
    seed=settings[:mcmcseed])

  # Dump output
  BSON.bson("$(results_dir)/output.bson", out)
  BSON.bson("$(results_dir)/simdat.bson", Dict(:simdat => config[:simdat]))

  println("Completed!")
end  # simfn
