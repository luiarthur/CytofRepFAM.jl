println("Testing...")

using Test
using CytofRepFAM
using BSON
using Random
using Distributions

# tests for feature allocation model
include("fam_tests.jl")

# tests for repulsive feature allocation model with feature selection
include("repfamFS_tests.jl")
