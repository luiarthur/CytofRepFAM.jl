module PlotUtils

include(joinpath(@__DIR__, "imports.jl"))

# Load python defs
path_to_plot_defs = @__DIR__
pushfirst!(PyCall.PyVector(PyCall.pyimport("sys")."path"), "$path_to_plot_defs")
plot_yz = PyCall.pyimport("plot_yz")
blue2red = PyCall.pyimport("blue2red")
pyrange(n) = collect(range(0, stop=n-1))

# General plot settings
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 15
rcParams["xtick.labelsize"] = 15
rcParams["ytick.labelsize"] = 15
rcParams["figure.figsize"] = (6, 6)


function quantiles(X, q; dims, drop=false)
  Q = mapslices(x -> quantile(x, q), X, dims=dims)
  out = drop ? dropdims(Q, dims=dims) : Q
  return out
end

extract(s, out) = [o[s] for o in out]

function boxplot(x; showmeans=true, whis=[2.5, 97.5], showfliers=false, kw...)
  plt.boxplot(x, showmeans=showmeans, whis=whis, showfliers=showfliers; kw...)
end

function add_gridlines_Z(Z)
  J, K = size(Z)
  for j in pyrange(J)
    plt.axhline(y=j+.5, color="grey", linewidth=.5)
  end

  for k in pyrange(K)
    plt.axvline(x=k+.5, color="grey", linewidth=.5)
  end
end

axhlines(x; kw...) = for xi in x plt.axhline(xi; kw...) end
cm_greys = plt.cm.get_cmap("Greys", 5)

# TODO
# Steal things from: Cytof5/sims/repfam_fs/test

function plot_Z(Z; colorbar=true)
  J, K = size(Z)
  p = plt.imshow(Z, aspect="auto", vmin=0, vmax=1, cmap=cm_greys)
  add_gridlines_Z(Z)
  plt.yticks(pyrange(J), pyrange(J) .+ 1, fontsize=rcParams["font.size"])
  plt.xticks(pyrange(K), pyrange(K) .+ 1, fontsize=rcParams["font.size"],
             rotation=90)
  if colorbar
    plt.colorbar()
  end
  return p
end

# plot y/z
function make_yz(y, Zs, Ws, lams, imgdir; vlim, 
                 w_thresh=.01, lw=3,
                 Z_true=nothing, 
                 markernames=[],
                 fs_y=rcParams["font.size"],
                 fs_z=rcParams["font.size"],
                 fs_ycbar=rcParams["font.size"],
                 fs_zcbar=rcParams["font.size"])
  # Make img dir if needed
  mkpath(imgdir)

  I = length(y)
  for i in 1:I
    idx_best = estimate_ZWi_index(Zs, Ws, i)

    Zi = Int.(Zs[idx_best])
    Wi = Float64.(Ws[idx_best][i, :])
    lami = Int64.(lams[idx_best][i])

    yi = Float64.(y[i])

    # plot Yi, lami
    plt.figure(figsize=(6, 6))
    plot_yz.plot_y(yi, Wi, lami, vlim=vlim, cm=blue2red.cm(9), lw=lw,
                   fs_xlab=fs_y, fs_ylab=fs_y, fs_lab=fs_y, fs_cbar=fs_ycbar,
                   markernames=markernames)
    plt.savefig("$(imgdir)/y$(i).pdf", bbox_inches="tight")
    plt.close()

    # plot Zi, Wi
    plt.figure(figsize=(6, 6))
    plot_yz.plot_Z(Zi, Wi, lami, w_thresh=w_thresh, add_colorbar=false,
                   fs_lab=fs_z, fs_celltypes=fs_z, fs_markers=fs_z,
                   fs_cbar=fs_zcbar, markernames=markernames)
    plt.savefig("$(imgdir)/Z$(i).pdf", bbox_inches="tight")
    plt.close()
  end

  if Z_true != nothing
    # plot Z true
    plt.figure(figsize=(6, 6))
    plot_yz.plot_Z_only(Z_true, fs=fs_z,
                        xlab="cell subpopulations", ylab="markers",
                        rotate_xticks=false)
    plt.savefig("$(imgdir)/Z_true.pdf", bbox_inches="tight")
    plt.close()

    # plot ZT true
    plt.figure(figsize=(6, 6))
    plot_yz.plot_Z_only(Z_true', fs=fs_z,
                        xlab="markers", ylab="cell subpopulations")
    plt.savefig("$(imgdir)/ZT_true.pdf", bbox_inches="tight")
    plt.close()
  end
end

grepKmcmc(s) = match(r"(?<=Kmcmc_)\d+", s)

function make_metrics(different_K_runs_dir, outputfname; thresh,
                      grepKmcmcfn=grepKmcmc)
  metrics = Dict{Int, Any}()

  for (root, dirs, files) in walkdir(different_K_runs_dir)
    for file in files
      if file == outputfname
        # Get path to output
        path_to_output = joinpath(root, file)
        println(path_to_output)

        # Parse Kmcmc
        Kmcmc = parse(Int, grepKmcmcfn(splitdir(root)[2]).match)

        # Load output
        output = BSON.load(path_to_output)

        # Append metrics to metrics dictionary
        # metrics[Kmcmc] = output[:metrics]
        metrics[Kmcmc] = Dict{Symbol, Any}()
        metrics[Kmcmc][:LPML] = output[:LPML]
        metrics[Kmcmc][:DIC] = output[:DIC]

        # Posterior samples of W
        Ws = extract(:theta__W, output[:samples][1])
 
        # Calibration metric
        cmetric = let
          Wmean = mean(Ws)
          metrics[Kmcmc][:I] = size(Wmean, 1)
          num_wik_lt_thresh = sum(Wmean .< thresh)
        end

        # Append calibration metric to metrics dictionary
        metrics[Kmcmc][:cmetric] = cmetric
       
        # Posterior samples of R (I x K) -- number of active features / sample
        # Rs (I x num_mcmc_samples)
        Rs = hcat([vec(sum(W .> 0, dims=2)) for W in Ws]...)

        # Mean of R, by sample (vector of length I)
        Rmean = vec(mean(Rs, dims=2))

        # Mean of R, by sample (vector of length I)
        R_025 = quantiles(Rs, .025, dims=2, drop=true)

        # Mean of R, by sample (vector of length I)
        R_975 = quantiles(Rs, .975, dims=2, drop=true)

        # Append Rmean to metrics dictionary
        metrics[Kmcmc][:Rmean] = Rmean

        # Append Rmean to metrics dictionary
        metrics[Kmcmc][:R_025] = R_025

        # Append Rmean to metrics dictionary
        metrics[Kmcmc][:R_975] = R_975
      end
    end
  end

  # All Kmcmcs
  Kmcmcs = sort(collect(keys(metrics)))

  # Create dir for metrics output
  metrics_dir = joinpath(different_K_runs_dir, "metrics")
  mkpath(metrics_dir)

  # Plot Calibration metric
  plt.plot([metrics[K][:cmetric] for K in Kmcmcs],
           [metrics[K][:LPML] for K in Kmcmcs], marker="o")
  for K in Kmcmcs
    plt.text(metrics[K][:cmetric], metrics[K][:LPML], K)
  end
  plt.xlabel(L"number of $W_{i,k}$ < 10%")
  plt.ylabel("LPML")
  plt.savefig(joinpath(metrics_dir, "calibration_metric.pdf"),
              bbox_inches="tight")
  plt.close()

  # Plot LPML, DIC
  for metric in [:LPML, :DIC]
    plt.plot(Kmcmcs, map(K -> metrics[K][metric], Kmcmcs), marker="o")
    plt.xlabel("K")
    plt.ylabel(String(metric))
    plt.savefig(joinpath(metrics_dir, "$(metric).pdf"), bbox_inches="tight")
    plt.close()
  end

  # Plot R
  I = metrics[Kmcmcs[1]][:I]
  plt.figure()
  for i in 1:I
    plt.subplot(I, 1, i)
    plt.plot(Kmcmcs, map(K -> metrics[K][:Rmean][i], Kmcmcs), marker="o")
    plt.fill_between(Kmcmcs,
                     map(K -> metrics[K][:R_025][i], Kmcmcs),
                     map(K -> metrics[K][:R_975][i], Kmcmcs), alpha=.5)
    if i == I
      plt.xlabel("K")
    end
    plt.ylabel("R$(i)")
  end
  plt.savefig(joinpath(metrics_dir, "R.pdf"), bbox_inches="tight")
  plt.close()

  return metrics
end

end  # PlotUtils
