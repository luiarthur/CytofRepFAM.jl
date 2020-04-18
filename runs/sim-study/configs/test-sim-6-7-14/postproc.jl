include("../../../PlotUtils/PlotUtils.jl")
include("../../../PlotUtils/imports.jl")

using Distributed
rmprocs(filter(w -> w > 1, workers()))
addprocs(27)

if length(ARGS) == 0
  results_dir = "/scratchdata/alui2/cytof/results/repfam/test-sim-6-7-14/"
else
  results_dir = ARGS[1]  # path to results directory
end

@everywhere include("../../../PlotUtils/PlotUtils.jl")
@everywhere include("../../../PlotUtils/imports.jl")

@everywhere function plot_params(samples, simdat, imgdir)
  # Create a directory for images / txt if needed.
  mkpath("$(imgdir)/txt/")

  # Extract samples
  extract(s) = map(o -> o[s], samples)
  Zs = extract(:theta__Z)
  Ws = extract(:theta__W)
  lams = extract(:theta__lam)
  alphas = extract(:theta__alpha)
  vs = extract(:theta__v)
  deltas = extract(:theta__delta)
  sig2s = extract(:theta__sig2)
  etas = extract(:theta__eta)
  omegas = extract(:omega)
  ps = [MCMC.sigmoid.(o) for o in omegas]

  # Number of samples
  I = length(sig2s[1])

  # Print Zmean
  open("$(imgdir)/txt/Zmean.txt", "w") do io
    writedlm(io, mean(Zs))
  end

  # Print Wmean
  open("$(imgdir)/txt/Wmean.txt", "w") do io
    writedlm(io, mean(Ws))
  end

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

  # Plot parameters
  PlotUtils.plot_W(Ws, imgdir=imgdir, W_true=simdat[:W])
  PlotUtils.plot_v(vs, imgdir)
  PlotUtils.plot_alpha(alphas, imgdir)
  PlotUtils.plot_mus(deltas, imgdir)
  PlotUtils.plot_sig2(sig2s, imgdir, sig2_true=simdat[:sig2])
  PlotUtils.plot_p(ps, imgdir)

  # Plot y / Z
  PlotUtils.make_yz(simdat[:y], Zs, Ws, lams, imgdir, vlim=(-4,4),
                    Z_true=simdat[:Z], w_thresh=0.0)
end

@everywhere function makeplots(path_to_output)
  output_dir = splitdir(path_to_output)[1]
  println(path_to_output)

  # Create img directory
  imgdir = joinpath(output_dir, "img")
  mkpath("$(imgdir)/txt")

  # Load simulated data
  path_to_simdat = joinpath(output_dir, "simdat.bson")
  simdat = BSON.load(path_to_simdat)[:simdat]

  # Load outputs
  output = BSON.load(path_to_output)
  samples = output[:samples][1]

  # Plot parameters
  plot_params(samples, simdat, imgdir)

  # Number of burn in iterations
  nburn = output[:nburn]

  # Plot loglikelihood
  println("Loglikes ...")
  if :ll in keys(output)
    ll = output[:ll]
    PlotUtils.plot_loglike(ll, imgdir, fname="loglike.pdf")
    PlotUtils.plot_loglike(ll[(nburn + 1):end],
                           imgdir, fname="loglike_postburn.pdf")
  end

  # Extract samples
  println("Extracting samples ...")
  extract(s) = map(o -> o[s], samples)
  Zs = extract(:theta__Z)
  Ws = extract(:theta__W)
  deltas = extract(:theta__delta)
  sig2s = extract(:theta__sig2)
  etas = extract(:theta__eta)

  # Number of samples
  I = length(sig2s[1])

  # Plot missing mechanism
  println("Plotting missmechs ...")
  for i in 1:I
    PlotUtils.plot_missmech(output[:c].constants.beta, i, xlim=[-5, 0])
    plt.savefig("$(imgdir)/missmech_$(i).pdf", bbox_inches="tight")
    plt.close()
  end

  # Print missing mechanism
  open(joinpath(imgdir, "txt/beta.txt"), "w") do io
    writedlm(io, output[:c].constants.beta, ',')
  end

  # Plot data density
  println("Plotting dden ...")
  PlotUtils.plot_dden(ddens=output[:dden],
                      etas=etas, Ws=Ws, Zs=Zs, sig2s=sig2s, deltas=deltas,
                      ygrid=output[:c].constants.y_grid,
                      imgdir=imgdir, simdat=simdat)
end

### MAIN ###

# Name of output file
OUTPUT_FILE = "output.bson"

# Path to all output files
output_paths = [joinpath(root, OUTPUT_FILE)
                for (root, _, files) in walkdir(results_dir)
                if OUTPUT_FILE in files]

# Make plots in parallel
error_msg = pmap(makeplots, output_paths, on_error=identity)
println(error_msg)

println("DONE!")
