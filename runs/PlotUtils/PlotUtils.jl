module PlotUtils

include(joinpath(@__DIR__, "imports.jl"))
include(joinpath(@__DIR__, "dden_complete.jl"))

# Load python defs
path_to_plot_defs = @__DIR__
pushfirst!(PyCall.PyVector(PyCall.pyimport("sys")."path"), "$path_to_plot_defs")
plot_yz = PyCall.pyimport("plot_yz")
blue2red = PyCall.pyimport("blue2red")
Population = PyCall.pyimport("Population").Population
pyrange(n) = collect(range(0, stop=n-1))

function quantiles(X, q; dims, drop=false)
  Q = mapslices(x -> quantile(x, q), X, dims=dims)
  out = drop ? dropdims(Q, dims=dims) : Q
  return out
end

skipnan(x) = x[.!isnan.(x)]

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

function plot_missmech(beta, i; xlim=[-10, 5],
                      ygrid=range(xlim[1], stop=xlim[2], length=300))
  p = [CytofRepFAM.Model.prob_miss(y, beta[:, i]) for y in ygrid]
  plt.plot(ygrid, p, lw=3)
  plt.title("Missing Mechanism for Sample $(i)")
  plt.xlabel("Probability of Missing")
  plt.ylabel("Expression Level")
end

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

  # Make img/txt dir
  txtdir = "$(imgdir)/txt" 
  mkpath(txtdir)

  I = length(y)
  population = Population()
  for i in 1:I
    idx_best = estimate_ZWi_index(Zs, Ws, i)

    Zi = Int.(Zs[idx_best])
    Wi = Float64.(Ws[idx_best][i, :])
    lami = Int64.(lams[idx_best][i])

    open("$(txtdir)/Z$(i)_best.txt", "w") do io
      writedlm(io, Zi)
    end

    open("$(txtdir)/W$(i)_best.txt", "w") do io
      writedlm(io, Wi)
    end

    open("$(txtdir)/lam$(i)_best.txt", "w") do io
      writedlm(io, lami)
    end

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
                   fs_cbar=fs_zcbar, markernames=markernames,
                   population=population)
    plt.savefig("$(imgdir)/Z$(i).pdf", bbox_inches="tight")
    plt.close()

    # plot Zmean
    plot_yz.plot_Z_only(mean(Zs), fs=fs_z,
                        xlab="cell subpopulations", ylab="markers",
                        rotate_xticks=true)
    plt.colorbar()
    plt.savefig("$(imgdir)/Zmean.pdf", bbox_inches="tight")
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

function plot_loglike(loglike, imgdir; fname="loglike.pdf", title=nothing)
  plt.plot(loglike)
  plt.xlabel("iter")
  plt.ylabel("log-likelihood")
  if title != nothing
    plt.title(title)
  end
  plt.savefig("$(imgdir)/$(fname)", bbox_inches="tight")
  plt.close()
end

function plot_alpha(alphas, imgdir; printmean=true)
  plt.hist(alphas, density=true)
  plt.xlabel("alpha")
  plt.ylabel("density")
  plt.savefig("$(imgdir)/alpha.pdf", bbox_inches="tight")
  plt.close()

  # Print mean alpha
  if printmean
    open("$(imgdir)/txt/alpha_mean.txt", "w") do io
      writedlm(io, mean(alphas))
    end
  end
end

function plot_v(vs, imgdir; printmean=true)
  v_mat = hcat(vs...)

  if printmean
    open("$(imgdir)/txt/v_mean.txt", "w") do io
      writedlm(io, mean(vs))
    end
  end

  boxplot(v_mat')
  plt.xlabel("cell phenotypes")
  plt.ylabel("v")
  plt.savefig("$(imgdir)/v.pdf", bbox_inches="tight")
  plt.close()
end

function plot_p(ps, imgdir; printmean=true)
  ps_mat = Matrix(hcat(ps...)')

  # Plot boxplots
  plt.boxplot(ps_mat)
  plt.xlabel("sample")
  plt.ylabel("p")
  plt.savefig("$(imgdir)/p.pdf", bbox_inches="tight")
  plt.close()

  # Traceplot of p
  pcols = size(ps_mat, 2)
  for i in 1:pcols
    plt.plot(ps_mat[:, i])
    plt.savefig("$(imgdir)/p$(i)_trace.pdf", bbox_inches="tight")
    plt.close()
  end

  if printmean
    open("$(imgdir)/txt/p_mean.txt", "w") do io
      writedlm(io, mean(ps))
    end
  end
end

function plot_mus(deltas, imgdir; printmean=true)
  mus0 = Matrix(hcat([-cumsum(d[0]) for d in deltas]...)')
  mus1 = Matrix(hcat([cumsum(d[1]) for d in deltas]...)')

  if printmean
    open("$(imgdir)/txt/mus0_mean.txt", "w") do io
      writedlm(io, mean(mus0, dims=1))
    end

    open("$(imgdir)/txt/mus1_mean.txt", "w") do io
      writedlm(io, mean(mus1, dims=1))
    end
  end

  mus = [mus0 mus1]
  boxplot(mus)
  plt.axhline(0)
  plt.savefig("$(imgdir)/mus.pdf", bbox_inches="tight")
  plt.close()

  # Traceplot of mus
  plt.plot(mus0)
  plt.plot(mus1)
  plt.savefig("$(imgdir)/mus_trace.pdf", bbox_inches="tight")
  plt.close()
end

function plot_sig2(sig2s, imgdir; sig2_true=nothing, printmean=true)
  sig2s = Matrix(hcat(sig2s...)')

  if printmean
    open("$(imgdir)/txt/sig2_mean.txt", "w") do io
      writedlm(io, mean(sig2s, dims=1))
    end
  end

  boxplot(sig2s)
  if sig2_true != nothing
    axhlines(sig2_true)
  end
  plt.savefig("$(imgdir)/sig2.pdf", bbox_inches="tight")
  plt.close()

  # Traceplot of sig2
  plt.plot(sig2s)
  plt.savefig("$(imgdir)/sig2_trace.pdf", bbox_inches="tight")
  plt.close()
end

function plot_W(Ws, imgdir; W_true=nothing, xlabel="cell subpopulations")
  plot_W(Ws, imgdir=imgdir, W_true=W_true, xlabel=xlabel)
end

function plot_W(Ws; imgdir, W_true=nothing, xlabel="cell subpopulations")
  I, K = size(Ws[1])
  # NOTE: Julia bug. Splat (...) throws errors if number of elems is too large.
  # Workaround: use list comprehensions instead.
  # See: https://github.com/JuliaLang/julia/issues/30796
  # Ws_mat = cat(Ws..., dims=3)
  S = length(Ws)

  plt.figure()
  for i in 1:I
    plt.subplot(I, 1, i)
    # boxplot(Ws_mat[i, :, :]')
    Wis = [Ws[s][i, k] for s in 1:S, k in 1:K]
    boxplot(Wis)
    plt.ylabel("W$(i)")
    if W_true != nothing
      K_true = size(W_true, 2)
      for k in 1:K_true
        plt.axhline(W_true[i, k])
      end
    end
  end
  plt.xlabel(xlabel)
  plt.savefig("$(imgdir)/W.pdf", bbox_inches="tight")
  plt.close()
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

function plot_dden(; ddens, etas, Ws, Zs, sig2s, deltas,
                   ygrid, imgdir, printmean=true, 
                   plot_true=true,
                   simdat=nothing, xlabel="expression level", ylabel="density",
                   dfs=nothing,
                   dden_xlim=(-6, 6))
  # Make directories if needed
  mkpath("$(imgdir)/dden")
  mkpath("$(imgdir)/txt")

  # Cache data sizes
  I, J = size(ddens[1])

  # split etas into etas0, etas1
  etas0 = [x[0] for x in etas]
  etas1 = [x[1] for x in etas]

  # Split mus from deltas
  mus0 = Matrix(hcat([-cumsum(d[0]) for d in deltas]...)')
  mus1 = Matrix(hcat([cumsum(d[1]) for d in deltas]...)')

  for i in 1:I
    for j in 1:J
      print("\r i: $i j: $j  ")

      # Get data density for sample i, marker j.
      # NOTE: `dden_ij` will be a vector (of length `length(ddens)`) of vectors
      # (of length `length(ygrid)`), and represents the posterior density at
      # the corresponding ygrid values.
      dden_ij = [ddij[i, j] for ddij in ddens]

      # 95% CI lower bound for dden ij
      dden_ij_lower = [quantile([ddij[g] for ddij in dden_ij], .025)
                       for g in 1:length(ygrid)]

      # 95% CI upper bound for dden ij
      dden_ij_upper = [quantile([ddij[g] for ddij in dden_ij], .975)
                       for g in 1:length(ygrid)]

      # Plot object of the credible interval for dden[i, j]
      p_ci_obs = plt.fill_between(ygrid, dden_ij_lower, dden_ij_upper,
                                  alpha=.5, color="blue")

      plt.xlabel(xlabel)
      plt.ylabel(ylabel)
      if dden_xlim != nothing
        plt.xlim(dden_xlim)
      end

      # Add eta to posterior dden plot
      eta0_ij_mean = mean(eta[0][i, j, :] for eta in etas)
      L0 = length(eta0_ij_mean)
      plt.scatter(mean(mus0, dims=1), zeros(L0), 
                  s=eta0_ij_mean * 60 .+ 10, marker="X", color="green")

      eta1_ij_mean = mean(eta[1][i, j, :] for eta in etas)
      L1 = length(eta1_ij_mean)
      p_mu = plt.scatter(mean(mus1, dims=1), zeros(L1),
                         s=eta1_ij_mean * 60 .+ 10, marker="X",
                         color="green")

      if printmean
        mkpath("$(imgdir)/txt/eta")

        open("$(imgdir)/txt/eta/eta0_i$(i)_j$(j)_mean.txt", "w") do io
          writedlm(io, eta0_ij_mean)
        end

        open("$(imgdir)/txt/eta/eta1_i$(i)_j$(j)_mean.txt", "w") do io
          writedlm(io, eta1_ij_mean)
        end
      end

      if simdat != nothing
        # Plot simulated data truth
        p_yobs = sns.kdeplot(skipnan(simdat[:y][i][:, j]),
                             color="red", bw=.1, label="tmp")
        p_yobs = p_yobs.get_legend_handles_labels()[1][1]

        # Histogram of observed data only
        # plt.hist(skipnan(simdat[:y][i][:, j]), color="red",
        #          alpha=.3, label="y (observed)",
        #          density=true, bins=30)

        if :eta in keys(simdat)
          eta_true = simdat[:eta]
        else
          eta_true = Dict(0 => ones(I, J, simdat[:L][0]),
                          1 => ones(I, J, simdat[:L][1]))
        end

        dgrid = dden_complete(ygrid, simdat[:W], eta_true,
                              simdat[:Z], simdat[:mus],
                              simdat[:sig2], i=i, j=j)

        if plot_true
          p_truth_complete, = plt.plot(ygrid, dgrid, color="grey", ls="--")
        end
      else
        # Plot histogram of observed data
        p_yobs = sns.kdeplot(skipnan(dfs.data.y[i][:, j]),
                             color="red", bw=.1, label="tmp")
        p_yobs = p_yobs.get_legend_handles_labels()[1][1]
      end

      # Number of MCMC samples
      num_mcmc_samples = length(Ws)

      # Create mus (vector)
      mus = [Dict(0 => -cumsum(delta[0]),
                  1 => cumsum(delta[1]))
             for delta in deltas]

      # Plot complete posterior density (obs and imputed)
      dd_complete_post = [dden_complete(ygrid,
                                        Ws[b],
                                        etas[b],
                                        Zs[b],
                                        mus[b],
                                        sig2s[b], i=i, j=j)
                          for b in 1:num_mcmc_samples]

      dd_complete_post = hcat(dd_complete_post...)  # legnth(ygird) x num_mcmc_samples
      dcp_lower = quantiles(dd_complete_post, .025, dims=2, drop=true)
      dcp_upper = quantiles(dd_complete_post, .975, dims=2, drop=true)
      p_ci_complete = plt.fill_between(ygrid, dcp_lower, dcp_upper, alpha=.5,
                                       color="orange")

      # TODO: PICK UP HERE
      if plot_true
        plt.legend([p_truth_complete, p_ci_complete,
                    p_yobs, p_ci_obs, p_mu],
                   ["truth (complete)", "95% CI (complete)",
                    "y (obs)", "95% CI (obs)", PyPlot.L"$\mu^\star$"])
      else
        plt.legend([p_ci_complete,
                    p_yobs, p_ci_obs, p_mu],
                   ["95% CI (complete)",
                    "y (obs)", "95% CI (obs)", PyPlot.L"$\mu^\star$"])
      end
      plt.savefig("$(imgdir)/dden/dden_i$(i)_j$(j).pdf",
                  bbox_inches="tight")
      plt.close()
    end
  end
  println()
end

end  # PlotUtils
