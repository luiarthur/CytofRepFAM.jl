import Dates
println(Dates.now())

include(joinpath(@__DIR__, "run_defs.jl"))

using Distributed
using Random
using Distributions
using BSON
using DelimitedFiles

# Simulation name.
simname = basename("$(@__DIR__)")
println(simname); flush(stdout)

# NOTE: write to scratchdir
function outdir_suffix(pmiss, phi, zind)
  return "pmiss$(pmiss)-phi$(phi)-zind$(zind)"
end

# read Commandline arguments.
RESULTS_DIR_PREFIX = ARGS[1]
AWS_BUCKET_PREFIX = ARGS[2]
PROPMISSINGSCALE = parse(Float64, ARGS[3])
PHI = parse(Float64, ARGS[4])
ZIND = parse(Int, ARGS[5])
ISTEST = parse(Int, ARGS[6]) > 0  # this argument should be 0 or 1 (for test)

# Create settings dictionary.
settings = let
  Nfac = 2000
  temperatures = [1.0, 1.003, 1.006, 1.01]
  Z = let
    path_to_Z = joinpath(@__DIR__, "Z$(ZIND).txt")
    readdlm(path_to_Z, Int, comments=true)
  end
  W = readdlm(joinpath(@__DIR__, "W.txt"), comments=true)
  _outdir_suffix = outdir_suffix(PROPMISSINGSCALE, Int(PHI), ZIND)
  results_dir = "$(RESULTS_DIR_PREFIX)/$(_outdir_suffix)"
  aws_bucket = "$(AWS_BUCKET_PREFIX)/$(_outdir_suffix)"
  if ISTEST
    nsamps = 10
    nburn = 10
  else
    # nsamps = 1000
    # nburn = 4000
    nsamps = 3000
    nburn = 6000
  end

  Dict(:Nfac => Nfac,
       :N => [1, 1] * Nfac,
       :thin_samps => 1,
       :Lmcmc => Dict(0 => 3, 1 => 3),
       :Kmcmc => 15,
       :pthin => 5,
       :dataseed => 1,
       :mcmcseed => 1,
       :temperatures => temperatures,
       :ntemps => length(temperatures),
       :Z => Z,
       :W => W,
       :results_dir => results_dir,
       :aws_bucket => aws_bucket,
       :nsamps => nsamps,
       :nburn => nburn,
       :propmissingscale => PROPMISSINGSCALE,
       :phi => PHI,
       :batchprop => 0.05,
       :dataseed => 1,
       :mcmcseed => 1)
end

# Set number of cores
ncores = settings[:ntemps]
addprocs(ncores)
println("nprocs: $(nprocs())")
println("nworkers: $(nworkers())")
println("system nproc: $(parse(Int, read(`nproc`, String)))")

# Import packages on all workers
@everywhere include("imports.jl")  # Load on all nodes.


println("Starting job ..."); flush(stdout)
println(Dates.now())

sim_and_run(settings)

println("DONE with job!"); flush(stdout)
println(Dates.now())
