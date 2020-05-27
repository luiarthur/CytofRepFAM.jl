ENV["PYCALL_JL_RUNTIME_PYTHON"] = "venv/bin/python"
import Pkg; Pkg.activate(joinpath(@__DIR__, "../../../../"))


using BSON
using PyCall
using PyPlot
using Pandas
using Seaborn
const pd = Pandas
const sns = Seaborn

const TSNE = pyimport("sklearn.manifold").TSNE

# Generate data
N = 200
K = 5
Y = vcat(randn(N, K), randn(N, K) .+ 4)
df = pd.DataFrame(Y, columns="V" .* string.(1:K))
df["class"] = vcat(fill(1, N), fill(2, N))

# Fit TSNE
tsne = TSNE(verbose=2, random_state=0)
@time tsne.fit(Y)

sns.pairplot(x_vars=:V1, y_vars=:V2, data=df, hue=:class,
             plot_kws=edgecolor=Dict(:linewidth => 0, :s=> 13),
             aspect=1, height=5)
plt.savefig("test.pdf", bbox_inches="tight")
plt.close()


path_to_simdat = ""
