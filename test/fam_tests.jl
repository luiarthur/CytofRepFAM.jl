#= Load these if running these tests in terminal
import Pkg; Pkg.activate("../")  # CytofRepFAM
using CytofRepFAM, Random, BSON, Test, Distributions
=#

printstyled("Test fitting FAM on simulated data ...\n", color=:yellow)
@testset "Fitting FAM on simulated data." begin
  # Simulation truth:
  @time dat = let
    J = 8  # number of markers
    N = [3, 1, 2] * 100  # number of observations
    I = length(N)  # number of samples
    K = 4  # number of latent features
    L = Dict(0 => 5, 1 => 3)  # number of density mixture components
    CytofRepFAM.Model.genData(J, N, J, L, sortLambda=true)
  end

  # Data to feed into model
  y = CytofRepFAM.Model.Data(dat[:y])

  # Model constants:
  c = let
    Kmcmc = 4
    Lmcmc = Dict(0 =>5, 1 => 3)
    CytofRepFAM.Model.defaultConstants(y, Kmcmc, Lmcmc,
                                       noisyDist=Normal(0, sqrt(10)),
                                       yBounds=[-6., -3., -2.])
  end

  # Print model constants
  CytofRepFAM.Model.printConstants(c)

  # Initialize model state
 
  # With Mclust
  @time init = CytofRepFAM.Model.smartInit(c, y)

  # Without Mclust
  # @time init = CytofRepFAM.Model.genInitialState(c, y)

  @time output = CytofRepFAM.Model.cytof_fit(
    init, c, y,
    nmcmc=20, nburn=20,
    computeLPML=true,
    computeDIC=true,
    computedden=true,
    joint_update_Z=false)

  # Create directory for output if needed.
  output_dir = "results/fam/"
  mkpath(output_dir)

  # Save model output 
  BSON.bson("$(output_dir)/output.bson", output)

  @test true
end
