include("../../../PlotUtils/PlotUtils.jl")
include("../../../PlotUtils/imports.jl")

if length(ARGS) == 0
  results_dir = "/scratchdata/alui2/cytof/results/repfam/test-sim-6-7-3"
else
  results_dir = ARGS[1]  # path to results directory
end

"""
Returns a list of monitors, where each element is the result from
one temperature.
"""
function getsamples(path_to_output)
  output = BSON.load(path_to_output)
  samples = output[:samples]
  num_tempers = length(output[:tempers])
  num_monitors = length(samples)

  return [let
     monitor = [let
                  num_samples = length(output[:samples][m])
                  [output[:samples][m][i][t] for i in 1:num_samples]
                end for m in 1:num_monitors]
   end for t in 1:num_tempers]
end

using Distributed
rmprocs(filter(w -> w > 1, workers()))
addprocs(20)

@everywhere include("../../../PlotUtils/PlotUtils.jl")
@everywhere include("../../../PlotUtils/imports.jl")

@everywhere function plot_params(samples, simdat, imgdir)
  # Create a directory for images / txt if needed.
  mkpath("$(imgdir)/txt/")

  # Extract samples
  extract(s) = map(o -> o[s], samples[1])
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

  # Plot y / Z
  PlotUtils.make_yz(simdat[:y], Zs, Ws, lams, imgdir, vlim=(-4,4),
                    Z_true=simdat[:Z], w_thresh=0.0)
end

function makeplots(path_to_output, imgdir)
  mkpath("$(imgdir)/txt")

  output_dir = splitdir(path_to_output)[1]
  path_to_simdat = joinpath(output_dir, "simdat.bson")
  simdat = BSON.load(path_to_simdat)[:simdat]

  println(path_to_output)
  output = BSON.load(path_to_output)
  samples = output[:samples][1]

  # Swap props
  println("Making swap proprs ...")
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
  println("Loglikes ...")
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

  # Extract samples
  println("Extracting samples ...")
  extract(s) = map(o -> o[s], samples[1])
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

# PATH TO ALL OUTPUT FILES
output_paths = [joinpath(root, OUTPUT_FILE)
                for (root, _, files) in walkdir(results_dir)
                if OUTPUT_FILE in files]

# Plot the posterior distributions of parameters
path_to_output = joinpath(results_dir, "maxtemp4-ntempts20-degree1-N500/output.bson")
path_to_simdat = joinpath(results_dir, "maxtemp4-ntempts20-degree1-N500/simdat.bson")
samples = getsamples(path_to_output);
simdat = BSON.load(path_to_simdat)[:simdat];
@everywhere simdat = $simdat
genimgdir(t) = joinpath(results_dir, "maxtemp4-ntempts20-degree1-N500/img/temper_$(t)")
@everywhere plot_params(s, imgdir) = plot_params(s, simdat, imgdir)
_ = pmap(plot_params, samples, genimgdir.(1:length(samples)))

# Make general plots
imgdir = joinpath(results_dir, "maxtemp4-ntempts20-degree1-N500/img")
makeplots(path_to_output, imgdir)

# Remove extra processors
# rmprocs(filter(w -> w > 1, workers()))

println("DONE!")

# Detailed Analysis
output = BSON.load(path_to_output);
last_states = output[:all_last_states];
tempers = output[:tempers]
num_tempers = length(tempers)

ll_fn(s, t) = CytofRepFAM.Model.compute_marg_loglike(s,
                                                     output[:c].constants,
                                                     output[:d].data, t)
i, j = 1, 20
ll_fn(last_states[i].theta, tempers[j]) -
ll_fn(last_states[i].theta, tempers[i])
ll_fn(last_states[j].theta, tempers[i]) -
ll_fn(last_states[j].theta, tempers[j])


temper_mat = [MCMC.WSPT.compute_log_accept_ratio(
  ll_fn, (last_states[i].theta, last_states[j].theta),
  (tempers[i], tempers[j]), verbose=3)
for i in 1:num_tempers, j in 1:num_tempers]
plt.imshow(exp.(temper_mat), vmin=.01, vmax=1.0)
plt.colorbar()
plt.savefig("tmp.pdf", bbox_inches="tight")
plt.close()
