include(joinpath(@__DIR__, "../PlotUtils/PlotUtils.jl"))
include(joinpath(@__DIR__, "../PlotUtils/imports.jl"))

result_dir = ARGS[1]  # path to results directory

### Create metrics ###

# Name of output file
OUTPUT_FILE = "output.bson"

println("Producing metrics ...")

# results_dirs = ["results/sim-runs-$(dsize)/mm0/"
#                 for dsize in ("small", "large")]
results_dirs = nothing # TODO

for d in results_dirs
  println(d)
  metrics = PlotUtils.make_metrics(d, OUTPUT_FILE, thresh=.01)
end


println("DONE!")
