include("../PlotUtils/PlotUtils.jl")
include("../PlotUtils/imports.jl")

# using Distributed
# rmprocs(filter(w -> w > 1, workers()))
# addprocs(2)


if length(ARGS) == 0
  results_dir = "/scratchdata/alui2/cytof/results/repfam/test-sim-6-7"
else
  results_dir = ARGS[1]  # path to results directory
end

# @everywhere include("../PlotUtils/PlotUtils.jl")
# @everywhere include("../PlotUtils/imports.jl")

# @everywhere function makeplots(path_to_output)
function makeplots(path_to_output)
  output_dir = splitdir(path_to_output)[1]
  path_to_simdat = joinpath(output_dir, "simdat.bson")
  simdat = BSON.load(path_to_simdat)[:simdat]

  println(path_to_output)
  output = BSON.load(path_to_output)

  # Extract samples
  extract(s) = map(o -> o[s], output[:samples][1])
  Zs = extract(:theta__Z)
  Ws = extract(:theta__W)
  lams = extract(:theta__lam)
  alphas = extract(:theta__alpha)
  vs = extract(:theta__v)
  deltas = extract(:theta__delta)
  sig2s = extract(:theta__sig2)
  etas = extract(:theta__eta)

  # Number of samples
  I = length(sig2s[1])

  # Create a directory for images / txt if needed.
  dir_to_output, _ = splitdir(path_to_output)
  imgdir = joinpath(dir_to_output, "img")
  mkpath("$(imgdir)/txt/")

  # Print R
  computeR(w) = vec(sum(w .> 0, dims=2))
  open("$(imgdir)/txt/Rcounts.txt", "w") do io
    I, K = size(Ws[1])
    Rs = map(computeR, Ws)
    Rcounts = zeros(Int, I, K)
    for R in Rs
      for i in 1:I
        Rcounts[i, R[i]] += 1
      end
    end
    writedlm(io, Rcounts)
  end

  # For parallel tempering
  if :swapcounts in keys(output) && :paircounts in keys(output)
    swapcounts = output[:swapcounts]
    paircounts = output[:paircounts]
    tempers = output[:tempers]
    ntempers = length(tempers)
    plt.imshow(swapcounts ./ (paircounts .+ 1e-6))
    plt.colorbar()
    plt.xticks(0:(ntempers-1), round.(tempers, digits=2), rotation=45)
    plt.yticks(0:(ntempers-1), round.(tempers, digits=2))
    plt.savefig("$(imgdir)/swapprops.pdf", bbox_inches="tight")
    plt.close()
  end

  # Number of burn in iterations
  nburn = output[:nburn]

  # Plot loglikelihood
  llkey = if :loglike in keys(output)
    PlotUtils.plot_loglike(output[:loglike], imgdir, fname="loglike.pdf")
    PlotUtils.plot_loglike(output[:loglike][(nburn + 1):end],
                           imgdir, fname="loglike_postburn.pdf")
  elseif :lls in keys(output)
    lls = output[:lls]
    for t in 1:ntempers
      temper = tempers[t]
      title = "temperature: $(temper)"
      ll_path = joinpath(imgdir, "loglike")
      mkpath(ll_path)
      PlotUtils.plot_loglike(lls[t], ll_path, fname="loglike_$(t).pdf",
                             title=title)
      PlotUtils.plot_loglike(lls[t][(nburn + 1):end],
                             ll_path, fname="loglike_postburn_$(t).pdf",
                             title=title)
    end
  end

  # Plot missing mechanism
  for i in 1:I
    PlotUtils.plot_missmech(output[:c].constants.beta, i, xlim=[-5, 0])
    plt.savefig("$(imgdir)/missmech_$(i).pdf", bbox_inches="tight")
    plt.close()
  end

  # Print missing mechanism
  open(joinpath(imgdir, "txt/beta.txt"), "w") do io
    writedlm(io, output[:c].constants.beta, ',')
  end

  # Plot parameters
  PlotUtils.plot_W(Ws, imgdir=imgdir, W_true=simdat[:W])
  PlotUtils.plot_v(vs, imgdir)
  PlotUtils.plot_alpha(alphas, imgdir)
  PlotUtils.plot_mus(deltas, imgdir)
  PlotUtils.plot_sig2(sig2s, imgdir, sig2_true=simdat[:sig2])

  # TODO: Plot data density
  PlotUtils.plot_dden(ddens=output[:dden],
                      etas=etas, Ws=Ws, Zs=Zs, sig2s=sig2s, deltas=deltas,
                      ygrid=output[:c].constants.y_grid,
                      imgdir=imgdir, simdat=simdat)

  # Plot y / Z
  PlotUtils.make_yz(simdat[:y], Zs, Ws, lams, imgdir, vlim=(-4,4),
                    Z_true=simdat[:Z], w_thresh=0.0)
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
status = map(makeplots, output_paths)
println(status)

# Remove extra processors
# rmprocs(filter(w -> w > 1, workers()))

println("DONE!")


