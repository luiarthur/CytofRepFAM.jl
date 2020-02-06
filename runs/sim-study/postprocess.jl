include("../PlotUtils/PlotUtils.jl")
include("../PlotUtils/imports.jl")

using Distributed
rmprocs(filter(w -> w > 1, workers()))
addprocs(10)


if length(ARGS) == 0
  results_dir = "/scratchdata/alui2/cytof/results/repfam/test-sim-6-5"
else
  results_dir = ARGS[1]  # path to results directory
end

@everywhere include("../PlotUtils/PlotUtils.jl")
@everywhere include("../PlotUtils/imports.jl")

@everywhere function makeplots(path_to_output)
  output_dir = splitdir(path_to_output)[1]
  path_to_simdat = joinpath(output_dir, "simdat.bson")
  simdat = BSON.load(path_to_simdat)[:simdat]

  println(path_to_output)
  output = BSON.load(path_to_output)
  extract(s) = map(o -> o[s], output[:samples][1])
  Zs = extract(:theta__Z)
  Ws = extract(:theta__W)
  lams = extract(:theta__lam)

  # Create a directory for images if needed.
  dir_to_output, _ = splitdir(path_to_output)
  imgdir = joinpath(dir_to_output, "img", )
  mkpath(imgdir)

  PlotUtils.make_yz(simdat[:y], Zs, Ws, lams, imgdir, vlim=(-4,4),
                    Z_true=simdat[:Z])
  PlotUtils.plotW(Ws, imgdir=imgdir, W_true=simdat[:W])
end

# # Get some simulation truths
# simdat_small = BSON.load("../data/simdat-nfac500.bson")
# simdat_large = BSON.load("../data/simdat-nfac5000.bson")
# Z_true_small = simdat_small[:Z]
# Z_true_large = simdat_large[:Z]
# W_true_small = simdat_small[:W]
# W_true_large = simdat_large[:W]
# 
# # Plot simulation truth Z (small data)
# PlotUtils.plot_yz.plot_Z_only(Z_true_small,
#                               fs=PlotUtils.rcParams["font.size"],
#                               xlab="cell subpopulations", ylab="markers")
# plt.savefig(joinpath(results_dir, "sim-runs-small/Z_true.pdf"))
# plt.close()
# 
# # Plot simulation truth Z (large data)
# PlotUtils.plot_yz.plot_Z_only(Z_true_large,
#                               fs=PlotUtils.rcParams["font.size"],
#                               xlab="cell subpopulations", ylab="markers")
# plt.savefig(joinpath(results_dir, "sim-runs-large/Z_true.pdf"))
# plt.close()
# 
# # Print simulation truth W (small data)
# open(joinpath(results_dir, "sim-runs-small/W_true.txt"), "w") do io
#   writedlm(io, W_true_small, ',')
# end
# 
# # Print simulation truth W (large data)
# open(joinpath(results_dir, "sim-runs-large/W_true.txt"), "w") do io
#   writedlm(io, W_true_large, ',')
# end
# 
### MAIN ###

# Name of output file
OUTPUT_FILE = "output.bson"

# PATH TO ALL OUTPUT FILES
output_paths = [joinpath(root, OUTPUT_FILE)
                for (root, _, files) in walkdir(results_dir)
                if OUTPUT_FILE in files]

# Make y, Z plots.
status = pmap(makeplots, output_paths)
println(status)

println("DONE!")
