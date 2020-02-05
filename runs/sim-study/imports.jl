println("Loading packages. This may take a few minutes ...")
flush(stdout)

# NOTE: These libs have to be imported in main process and worker processes.
import Pkg; Pkg.activate(joinpath(@__DIR__, "../../"))  # CytofRepFAM
using CytofRepFAM
using Random, Distributions, BSON

println("Finished loading packages. This may take a few minutes ...")
flush(stdout)
