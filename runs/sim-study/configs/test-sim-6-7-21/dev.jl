# Load imports
include(joinpath(@__DIR__, "../../imports.jl"))
include(joinpath(@__DIR__, "../../../PlotUtils/PlotUtils.jl"))
import StatsBase
const plt = PlotUtils.plt
import Seaborn
const sns = Seaborn

# Path to simulation with: Z2, MCMC2, phi1
path_to_results = "/scratchdata/alui2/cytof/results/repfam/test-sim-6-7-21/dataseed1-mcmcseed2-phi1.0-Zind2"

# Read data
output = BSON.load(joinpath(path_to_results, "output.bson"))

# Extract posterior
extract(s) = map(o -> o[s], output[:samples][1])

# Get some extracts
Zs = extract(:theta__Z)
Ws = extract(:theta__W)
lams = extract(:theta__lam)

# Get point estimate index
idx_best = [begin
              PlotUtils.estimate_ZWi_index(Zs, Ws, i)
            end for i in 1:output[:d].data.I]

# Best lambdas
lam_best = [lams[idx_best[i]][i] for i in 1:output[:d].data.I]

# Get counts of each feature in lam_1-est and lam_2-est
cm1 = StatsBase.countmap(lam_best[1])
cm2 = StatsBase.countmap(lam_best[2])

[k => v / output[:d].data.N[1] for (k, v) in cm1]
[k => v / output[:d].data.N[2] for (k, v) in cm2]


data = vcat([begin
  y = output[:lastState].theta.y_imputed[1][lam_best[1] .== k, 7:9]
  cluster = fill(k, size(y, 1))
  [y cluster]
end for k in 7:10]...)

df = DataFrame(data, [:m7, :m8, :m9, :feature])
df_stacked = stack(df, Not(:feature), variable_name=:marker, value_name=:y)
CSV.write("img/df.csv", df_stacked)

# Do the rest in python: dev.py
