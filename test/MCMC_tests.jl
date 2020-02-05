#= Load these if running these tests in terminal
import Pkg; Pkg.activate("../")  # CytofRepFAM
import CytofRepFAM.MCMC
=#

@testset "Test logit / sigmoid" begin
  @test MCMC.logit(.6) ≈ log(.6 / .4)
  @test MCMC.sigmoid(3.) ≈ 1 / (1 + exp(-3.))
end


@testset "Test logsumexp" begin
  x = randn(10)
  @test MCMC.logsumexp(x) ≈ log(sum(exp.(x)))
  @test MCMC.logsumexp(x[1:3]) ≈ MCMC.logsumexp(x[1], x[2], x[3])
end
