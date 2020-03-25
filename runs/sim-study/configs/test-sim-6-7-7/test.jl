include("../../imports.jl")  # Load on main node.

# NOTE: Change this.
SIMNAME = "$(basename(@__DIR__))-test"
SCRATCH_DIR = "/scratchdata/$(ENV["USER"])/cytof/results/repfam"
RESULTS_DIR_PREFIX = joinpath(SCRATCH_DIR, SIMNAME)
AWS_BUCKET_PREFIX = "s3://cytof-repfam/$(SIMNAME)"

# NOTE: read configs
include("settings.jl")

for setting in settings
  setting[:results_dir] = "$(RESULTS_DIR_PREFIX)/$(setting[:outdir_suffix])"
  setting[:aws_bucket] = "$(AWS_BUCKET_PREFIX)/$(setting[:outdir_suffix])"
end

println("Sourcing files ..."); flush(stdout)
include("simfn.jl")

function sim(setting)
  results_dir = setting[:results_dir]
  aws_bucket = setting[:aws_bucket]
  mkpath(results_dir)

  println("Running $(results_dir)")
  flush(stdout)

  CytofRepFAM.Model.redirect_stdout_to_file("$(results_dir)/log.txt") do
    simfn(setting)
  end

  # Send to S3.
  CytofRepFAM.Model.s3sync(from=results_dir, to=aws_bucket,
                           tags=`--exclude '*.nfs'`)

  return
end

# NOTE:
# - f::Function: A function which takes one argument of type Dict{Any}
# - settings::Vector{Dict}): A vector of settings
println("Starting jobs ..."); flush(stdout)
_ = sim(settings[1])
println("DONE with all runs!")
