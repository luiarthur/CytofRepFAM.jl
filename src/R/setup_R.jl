module Rplots

using RCall
const RSOURCE_DIR = "$(pwd())/src/R"
const RPLOT_SOURCE = "$(RSOURCE_DIR)/plots.R"
rplotlib() = R"source($RPLOT_SOURCE, chdir=TRUE)"

end # Rplots

#= TEST
using CytofRepFAM
CytofRepFAM.Rplots.RSOURCE_DIR
CytofRepFAM.Rplots.rplotlib()
=#
