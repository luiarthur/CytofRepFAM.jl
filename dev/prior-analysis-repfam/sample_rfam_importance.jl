import Pkg; Pkg.activate("../../")

using Distributions
using Statistics
using LinearAlgebra
using PyPlot
PyPlot.matplotlib.use("Agg")

function pairwise_dist(Z)
  K = size(Z, 2)
  return [norm(Z[:, k] - Z[:, h], 1) for k in 2:K for h in 1:(k - 1)]
end

Dmean(Z) = mean(pairwise_dist(Z))
Dmin(Z) = minimum(pairwise_dist(Z))

function rep_fn(Z, phi)
  d = pairwise_dist(Z)
  if phi == 0
    return 1.0
  else
    log_rep_term = phi * sum(log1p.(-exp.(-d)))
    rep_term = exp(log_rep_term)
    return rep_term
  end
end

"""
Approximate the average distance between columns of rFAM
E[D|phi] using importance sampling.
"""
function mean_D_rfam(N::Integer, J::Integer, v::Vector{Float64},
                     repulsive_fn::Function; D=Dmean, return_samples=false)
  K = length(v)
  Zs = [[v[k] > rand() for j in 1:J, k in 1:K] for _ in 1:N]

  weights = repulsive_fn.(Zs)
  numer = sum(D.(Zs) .* weights)
  denom = sum(weights)

  if return_samples
    return (Zs, weights)
  else
    return numer / denom  # if NaN is returned, perhaps more samples are needed.
  end
end

function mean_D_rfam(N::Integer, J::Integer, K::Integer, alpha::Float64,
                     repulsive_fn::Function; D=Dmean)
  Zs = [let
          v = rand(Beta(alpha / K), 1, K)
          v .> rand(J, K)
        end for _ in 1:N]

  weights = repulsive_fn.(Zs)
  numer = sum(D.(Zs) .* weights)
  denom = sum(weights)

  if denom == 0
    println("WARNING: weights are all 0!")
    return 0.0
  else
    return numer / denom
  end
end


# MAIN
phis = [0, 1, 10, 100]

# Dmean
@time dmean = [mean_D_rfam(100000, 15, fill(.5, 25), Z -> rep_fn(Z, phi), D=Dmean)
               for phi in phis];
foreach(z -> println("phi=$(z[1]) => $(z[2])"), zip(phis, dmean))

out = [mean_D_rfam(100000, 15, fill(.5, 25), Z -> rep_fn(Z, phi), D=Dmean,
                   return_samples=true) for phi in phis];

plt.figure()
Zss = [let
         Zs, w = out[i]
         wsample(Zs, w, 1000)
       end for i in 1:length(phis)]

mean_d = [mean(Dmean.(Zs)) for Zs in Zss]
min_d = [minimum(Dmean.(Zs)) for Zs in Zss]
max_d = [maximum(Dmean.(Zs)) for Zs in Zss]
q975_d = [quantile(Dmean.(Zs), 0.975) for Zs in Zss]
q025_d = [quantile(Dmean.(Zs), 0.025) for Zs in Zss]
          
plt.plot(phis, mean_d, lw=3, marker="o", label="mean")
plt.plot(phis, min_d, ls=":", marker="o", label="min")
plt.plot(phis, max_d, ls=":", marker="o", label="max")
plt.plot(phis, q025_d, ls="-.", marker="o", label="Q2.5")
plt.plot(phis, q975_d, ls="-.", marker="o", label="Q97.5")
plt.xticks(phis, phis)
plt.xlabel("phi")
plt.ylabel("Average pairwise column distance in Z")
plt.xscale("log")
plt.legend()
plt.savefig("bla.pdf", bbox_inches="tight")
plt.close()

# Dmin
@time dmin = [mean_D_rfam(100000, 15, fill(.5, 25), Z -> rep_fn(Z, phi), D=Dmin)
              for phi in phis];
foreach(z -> println("phi=$(z[1]) => $(z[2])"), zip(phis, dmin));

# D > 2
D_gt_two(Z) = all(pairwise_dist(Z) .>= 2)
@time dgt2 = [mean_D_rfam(100000, 15, fill(.5, 25), Z -> rep_fn(Z, phi), D=D_gt_two)
              for phi in phis];
foreach(z -> println("phi=$(z[1]) => $(z[2])"), zip(phis, dgt2));

# D > u
phis = collect(0:5:30)
J, K = (15, 25)
for u in 0:5
  D_ge_u(Z) = minimum(pairwise_dist(Z)) >= u
  @time dgeu = [mean_D_rfam(10000, J, fill(.5, K), Z -> rep_fn(Z, phi), D=D_ge_u)
                for phi in phis]
  foreach(z -> println("phi=$(z[1]) => $(z[2])"), zip(phis, dgeu));
  plt.plot(phis, dgeu, label="u=$(u)", marker="o")
end
plt.legend()
plt.title("J=$J, K=$K")
plt.ylabel("Pr(min(D(Z)) ≥ u | phi)")
plt.xlabel("phi")
plt.axhline(0.95, ls=":", c="grey")
plt.savefig("img/pu.pdf")
plt.close()
