#= Load these if running these tests in terminal
import Pkg; Pkg.activate("../")  # CytofRepFAM
using CytofRepFAM, Test, BSON, Random, Distributions
=#

@testset "Test MCMC.Metropolis" begin
  x = 1.0
  logprob(x) = logpdf(Normal(0, 1), x)
  jumpdist(x) = Normal(x, .5)
  logjump(y, x) = logpdf(jumpdist(x), y)
  randjump(x) = rand(jumpdist(x))
  metro = MCMC.Metropolis(x, 0.0)
  samps = Float64[]
  for i in 1:10000
    append!(samps,
            MCMC.update!(metro, logprob, logjump, randjump, verbose=0).state)
  end

  @test true
end


