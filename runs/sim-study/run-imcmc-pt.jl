import Dates
println(Dates.now())

using Distributed
include("imports.jl")  # Load on main node.


# NOTE: Change this.
if length(ARGS) == 4
  SIMDIR = ARGS[1]
  RESULTS_DIR_PREFIX = ARGS[2]
  AWS_BUCKET_PREFIX = ARGS[3]
  ISTEST = parse(Int, ARGS[4]) > 0  # this argument should be 0 or 1 (for test)
else
  println("Usage: julia parsim.jl <simdir> <results_dir_prefix> <aws_bucket_prefix> <istest>")
  exit(1)
end

function setnumcores(n::Int)
  children_procs = filter(w -> w > 1, workers())
  rmprocs(children_procs)
  addprocs(n)
end

# NOTE: read configs
include("$(SIMDIR)/settings.jl")

for setting in settings
  setting[:results_dir] = "$(RESULTS_DIR_PREFIX)/$(setting[:outdir_suffix])"
  setting[:aws_bucket] = "$(AWS_BUCKET_PREFIX)/$(setting[:outdir_suffix])"
  if ISTEST
    setting[:nsamps] = 10
    setting[:nburn] = 10
  end
end

ncores = 24
setnumcores(ncores)
@everywhere include("imports.jl")  # Load on main node.

println("Sourcing files ..."); flush(stdout)
@everywhere SIMDIR = $SIMDIR
@everywhere include("$(SIMDIR)/simfn.jl")

@everywhere function sim(setting)
  results_dir = setting[:results_dir]
  aws_bucket = setting[:aws_bucket]
  mkpath(results_dir)

  println("Running $(results_dir)")
  flush(stdout)

  # FIXME: Can't do `asyncmap` and redirect outputs.
  # CytofRepFAM.Model.redirect_stdout_to_file("$(results_dir)/log.txt") do
  #   simfn(setting)
  # end

  simfn(setting)

  # Send to S3.
  CytofRepFAM.Model.s3sync(from=results_dir, to=aws_bucket,
                           tags=`--exclude '*.nfs'`)

  # Remove results to save space.
  # rm(results_dir, recursive=true)
  return
end

# NOTE:
# - f::Function: A function which takes one argument of type Dict{Any}
# - settings::Vector{Dict}): A vector of settings
println("Starting jobs ..."); flush(stdout)
ntasks = 6  # 6 runs at a time, because each run does 4-core PT
@time _ = asyncmap(sim, settings, ntasks=ntasks);
println("DONE with all runs!")
println(Dates.now())
