include("../../../PlotUtils/PlotUtils.jl")
include(joinpath(@__DIR__, "../../imports.jl"))
using PyPlot
const plt = PlotUtils.plt

function tomu(delta)
  mu0 = -cumsum(delta[0])
  mu1 = cumsum(delta[1])
  return [mu0; mu1]
end

function etas_tomat(etas)
  eta0s = [eta[0] for eta in etas]
  eta1s = [eta[1] for eta in etas]
  return (cat(eta0s..., dims=4),
          cat(eta1s..., dims=4))
end

# Results directory
results_dir = "/scratchdata/alui2/cytof/results/repfam/test-sim-6-7-13/"

# Name of output file
OUTPUT_FILE = "output.bson"

# Path to all output files
output_paths = [joinpath(root, OUTPUT_FILE)
                for (root, _, files) in walkdir(results_dir)
                if OUTPUT_FILE in files]

# Load outputs
path_to_output = output_paths[2]
output = BSON.load(path_to_output)
samples = output[:samples][1]

extract(sym) = [s[sym] for s in samples]

# mu*
mus = [tomu(d) for d in extract(:theta__delta)]
mu_mat = Matrix(hcat(mus...)')

# sig2
sig2s = extract(:theta__sig2)
sig2_mat = Matrix(hcat(sig2s...)')

# W
Ws = extract(:theta__W)

# Z
Zs = extract(:theta__Z)

# R
Rs = hcat([vec(sum(W .> 0, dims=2)) for W in Ws]...)

# eta
eta0_mat, eta1_mat = etas_tomat(extract(:theta__eta))

# Plot eta vs mu
for i in 1:2, j in 1:20
  println((i, j))

  plt.subplot(2, 1, 1)
  plt.plot(mu_mat[:, 1], eta0_mat[i, j, 1, :], lw=0.5)
  plt.xlabel(L"$\mu^\star_{0, 1}$")
  plt.ylabel(latexstring("\\eta^0_{1, $(j), 1}"))

  plt.subplot(2, 1, 2)
  plt.plot(sig2_mat[:, i], eta0_mat[i, j, 1, :], lw=0.5)
  plt.xlabel(latexstring("\\sigma^2_{$(i)}"))
  plt.ylabel(latexstring("\\eta^0_{1, $(j), 1}"))

  plt.savefig("img/eta0_$(i)_$(j)_1-vs-mu.pdf", bbox_inches="tight")
  plt.close()
end

