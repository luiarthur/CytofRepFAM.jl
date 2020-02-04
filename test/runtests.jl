println("Testing...")

using Test
using CytofRepFAM
import CytofRepFAM.MCMC
using BSON
using Random
using Distributions

# tests for MCMC
include("MCMC_tests.jl")

# tests for feature allocation model (FAM)
include("fam_tests.jl")

# tests for repulsive FAM with feature selection
include("repfamFS_tests.jl")

# tests for variational inference with FAM
include("vb_fam_tests.jl")
