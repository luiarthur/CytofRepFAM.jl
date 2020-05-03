import Pkg; Pkg.activate(joinpath(@__DIR__, "../../../../../"))
using Distributions
using PyPlot

"""
Paper: Bayesian inference for finite mixtures of univariate and multivariate
       skew-normal and skew-t distributions, Biostatistics 2010.
skew (delta): a real number in (-1, 1)
"""
function rand_skewnormal(loc, scale, skew)
    z = rand(TruncatedNormal(0, 1, 0, Inf))
    return loc + scale * skew * z + scale * sqrt(1 - skew ^ 2) * randn()
end

function rand_skewnormal(loc, scale, skew, dims...)
    z = rand(TruncatedNormal(0, 1, 0, Inf), dims...)
    return loc .+ scale * skew * z + scale * sqrt(1 - skew ^ 2) * randn(dims...)
end

x = rand_skewnormal(1, .5, -.97, 100000)
plt.hist(x, bins=100);

mean(x .< 0)

mean(x), std(x)


