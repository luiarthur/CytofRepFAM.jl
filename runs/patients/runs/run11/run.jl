# NOTE: Carefully review all the NOTES here.
using Distributed

# This is needed to compile once before loading on all cores.
include(joinpath(@__DIR__, "imports.jl"))

# NOTE: We are using 4 temperatures, so using 4 cores.
addprocs(4)

@everywhere include(joinpath(@__DIR__, "imports.jl"))
@everywhere const rfam = CytofRepFAM.Model

@everywhere function preprocess(y, lower, upper)
  @assert lower < upper
  idx_keep = .!vec(any(y .< lower, dims=2) .|
                   any(y .> upper, dims=2))
  return y[idx_keep, :]
end

@everywhere function run(phi, path_to_data::Vector{String},
                         results_dir, aws_bucket;
                         nburn=10000, nsamps=5000, thin=1,
                         K=25, L=Dict(0=>6, 1=>3),
                         tempers=[1, 1.003, 1.006, 1.01],
                         batchsizes=[200, 200], pthin=5,
                         y_lower=-7, y_upper=4)

  @assert length(path_to_data) == length(batchsizes)

  # Print setup
  println("Start time: ", Dates.now())
  println("pid: $(getpid())")
  println("Threads: $(Threads.nthreads())")
  println("pwd: $(pwd())")
  println("Data used:")
  foreach(println, "  " .* path_to_data)
  flush(stdout)

  # TODO: Make y
  function read_and_format(path::String)
    df = CSV.read(path, missingstring="NA")
    mat = Matrix(coalesce.(df, NaN))
    return preprocess(mat, y_lower, y_upper), names(df)
  end
  data = read_and_format.(path_to_data)
  y = [d[1] for d in data]
  println("data sizes: ", size.(y))
  markers = let
    _markers = [d[2] for d in data]
    @assert length(unique(_markers)) == 1
    _markers[1]
  end
  println("markers: $(String.(markers))")
  
  # Remove markers
  num_markers = length(markers)
  marker_index_to_remove = [2, 3, 4, 6, 8, 9, 10, 11, 12, 18, 19, 20, 21, 23,
                            26, 30, 31]
  marker_index_to_keep = filter(j -> !(j in marker_index_to_remove),
                                1:num_markers)
  println("Good markers: $(String.(markers[marker_index_to_keep]))")
  println("Bad markers: $(String.(markers[marker_index_to_remove]))")
  y = [yi[:, marker_index_to_keep] for yi in y]
  println("Reduced data sizes: ", size.(y))

  # Results dir
  mkpath(results_dir)

  # number of temperatures
  ntemps = length(tempers)

  # Repulsive FAM penalty scale
  # NOTE: Repulsive function.
  log_repulsive_fn = rfam.log_repulsive_fn_generator(phi)
  use_repulsive = phi > 0

  function init_state_const_data(y; K, L)
    deltaz_prior = TruncatedNormal(1.0, 0.1, 0.0, Inf)
    d = rfam.Data(y)
    c = rfam.defaultConstants(d, K, L,
                              sig2_prior=InverseGamma(3, 1),
                              delta0_prior=deltaz_prior,
                              delta1_prior=deltaz_prior,
                              alpha_prior=Gamma(1, 1),  # shape, scale
                              # NOTE:  These 3 items are important
                              yQuantiles=[.0, .25, .5], 
                              pBounds=[.05, .8, .05],
                              y_grid=collect(range(-10,
                                                   stop=4,
                                                   length=100)),
                              probFlip_Z=1.0,
                              log_repulsive_fn=log_repulsive_fn)

    # NOTE: To spread weights evenly.
    c.eta_prior = Dict(z => Dirichlet(L[z], 1) for z in 0:1)

    # Initialize state.
    # NOTE: Kmeans init
    # states = [rfam.smartInit(c, d, modelNames="kmeans",
    #                          seed=1 + s, iterMax=100)
    #           for s in 1:ntemps]
    # NOTE: Mclust init
    # VII: spherical, unequal volume
    states = [rfam.smartInit(c, d, modelNames="VII", seed=1 + s,
                             allowZDupCols=false) for s in 1:ntemps]

    # Create tunenrs
    t = rfam.Tuners(d.y, c.K)

    # Create covriates matrix
    X = rfam.eye(Float64, d.I)

    # Create model constants
    cfs = rfam.ConstantsFS(c)

    # NOTE: Use informative priors for W* and omega, to control feature
    # selection.
    # cfs.W_star_prior = Gamma(1000, 0.1) # shape, scale. Works, but extreme?
    # cfs.W_star_prior = Gamma(100, 0.1) # shape, scale. Not working.
    # cfs.W_star_prior = Gamma(1/c.K, 1) # shape, scale. Not working.
    # cfs.W_star_prior = Gamma(1/c.K, c.K * 10) # shape, scale. Not working.
    # cfs.W_star_prior = Gamma(1, 100) # shape, scale. Not super helpful?
    # cfs.W_star_prior = Gamma(10, 1) # shape, scale. Try?
    # cfs.W_star_prior = Gamma(5, 1) # shape, scale. Try?
    # cfs.omega_prior = Normal(-5, 1) # similar to p ~ Beta(1, 99)

    # From run1
    cfs.omega_prior = Normal(-3, 1.3)  # similar to Beta(1, 9)
    cfs.W_star_prior = Gamma(1.0, 2) # shape, scale

    dfs = rfam.DataFS(d, X)
    sfss = [let
              sfs = rfam.StateFS{Float64}(s, dfs, omega_prior=cfs.omega_prior,
                                          W_star_prior=cfs.W_star_prior,
                                          eps_r=0.05)
              # Fix omega
              K0 = 5
              sfs.omega .= CytofRepFAM.MCMC.logit(K0/K)
              println("Initial omega: $(sfs.omega)")
              sfs
            end for s in states]
    tfs = rfam.TunersFS(t, states[1], X)

    return Dict(:dfs => dfs, :cfs => cfs, :sfss => sfss, :tfs => tfs, :y => y)
  end

  # Parameters to monitor
  monitor1 = [:theta__Z, :theta__v, :theta__alpha,
              :r, :theta__lam, :W_star, :theta__eta,
              :theta__W, :theta__delta, :theta__sig2]  # removed omega
  monitor2 = [:theta__y_imputed, :theta__gam]

  # MCMC Specs
  nsamps_to_thin(nsamps::Int, nmcmc::Int) = max(1, div(nmcmc, nsamps))
  thin_samps = thin  # Factor to thin the
                     # primary parameters
  mcmc_iter = nsamps * thin_samps  # Number of MCMC iterations

  # LPML / DIC are computed based on `mcmc_iter` samples

  # Configurations: priors, initial state, data, etc.
  config = init_state_const_data(y, K=K, L=L)

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
    fix=[:omega],
    batchsizes=batchsizes,
    prior_thin=pthin,
    tempers=tempers,
    randpair=1.0, swap_freq=1.0,
    Z_marg_lamgam=1.0,
    Z_marg_lamgam_decay_rate=100.0,
    Z_marg_lamgam_min=1.0,
    use_repulsive=use_repulsive,
    printFreq=1,
    seed=1,
    sb_ibp=false,
    computedden=true,
    computeDIC=true,
    computeLPML=true,
    thins=[thin, nsamps_to_thin(10, mcmc_iter)],
    ndden_samps=200,
    time_updates=false,
    save_all_states=true,
    verbose=4)

  # Dump output
  BSON.bson("$(results_dir)/output.bson", out)

  println("Completed!")
  println("End time: ", Dates.now())
  flush(stdout)

  # Send results to aws
  rfam.s3sync(from=results_dir, to=aws_bucket, tags=`--exclude '*.nfs'`)

  return
end  # simfn

### MAIN ###

# READ COMMAND ARGS
phi = parse(Float64, ARGS[1])
path_to_data = String.(split(ARGS[2], ","))
results_dir = ARGS[3]
aws_bucket = ARGS[4]
istest = Bool(parse(Int, ARGS[5]))

println("phi: $phi")
println("data path: $path_to_data")
println("results dir: $results_dir")
println("AWS bucket: $aws_bucket")
println("is test: $istest")
flush(stdout)

if istest
  nsamps = 20
  nburn = 10
else
  nsamps = 5000
  nburn = 10000
end

@time run(phi, path_to_data, results_dir, aws_bucket,
          nsamps=nsamps, nburn=nburn)
