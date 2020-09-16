import Pkg; Pkg.activate("../../")

import CytofRepFAM.Model.ImportanceSampling: exp_f_rfam, min_pairwise_dist_ge_u, rep_fn
using PyPlot
PyPlot.matplotlib.use("Agg")

function get_exp_f(phis, us; N, J, K)
  [[let
      println("phi: $phi, u: $u")
      exp_f_rfam(N=N, J=J, K=K, 
                 repulsive_fn=Z->rep_fn(Z, phi),
                 f=Z->min_pairwise_dist_ge_u(Z, u))
    end for u in us] for phi in phis]
end

function plot_results(phis, us, exp_f, savepath)
  for z in zip(phis, exp_f)
    phi, f = z
    plt.plot(us, f, label="ϕ = $phi", marker="o")
  end
  plt.xticks(fontsize=15)
  plt.yticks(fontsize=15)
  plt.xlabel("threshold " * L"d̲", fontsize=18)
  plt.ylabel(L"\Pr\left(\min_{1 ≤ k₁ < k₂ ≤ K}~d(\mathbf{z}_{k_1}, \mathbf{z}_{k_2}) ≥ d̲ \mid \phi \right)",
             fontsize=18)
  plt.axhline(0.95, ls=":", color="grey")
  plt.legend(loc="lower left")
  plt.savefig(savepath, bbox_inches="tight")
  plt.close()
end

# phis = 0:5:50
phis = 0:5:30
us = 0:5

exp_f_data = get_exp_f(phis, us, N=10000, J=15, K=25)
exp_f_sim = get_exp_f(phis, us, N=10000, J=21, K=15)

plt.figure(figsize=(6, 6))
plot_results(phis, us, exp_f_data, "img/data_prior_analysis.pdf")
plt.figure(figsize=(6, 6))
plot_results(phis, us, exp_f_sim, "img/sim_prior_analysis.pdf")
