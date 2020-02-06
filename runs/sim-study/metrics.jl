include(joinpath(@__DIR__, "../PlotUtils/PlotUtils.jl"))
include(joinpath(@__DIR__, "../PlotUtils/imports.jl"))

if length(ARGS) == 0
  result_dir = ARGS[1]  # path to results directory
else
  result_dir = "/scratchdata/alui2/cytof/results/repfam/test-sim-6-5"
end

### Create metrics ###

# Name of output file
OUTPUT_FILE = "output.bson"

println("Producing metrics ...")

# directories to loop through
results_dirs = ["$(result_dir)/seed_$(seed)/scale_$(scale)"
                for seed in 1:3
                for scale in (0, 1, 10)]

for d in results_dirs
  println(d)
  metrics = PlotUtils.make_metrics(d, OUTPUT_FILE, thresh=.01)
end

println("DONE!")
