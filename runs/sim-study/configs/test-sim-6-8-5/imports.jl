println("Loading packages. This may take a few minutes ...")
flush(stdout)

# NOTE: These libs have to be imported in main process and worker processes.
import Pkg; Pkg.activate(joinpath(@__DIR__, "../../../../"))  # CytofRepFAM
using CytofRepFAM
const rfam = CytofRepFAM.Model

println("Finished loading packages.")
flush(stdout)
