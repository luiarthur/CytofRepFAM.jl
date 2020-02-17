module ParTemp

import Pkg
Pkg.activate(joinpath(@__DIR__, "../../"))  # CytofRepFam
import CytofRepFAM.MCMC
using Distributions

function metropolisMV(curr, logprob, S)
  cand = rand(MvNormal(curr, S))
  logU = log(rand())
  logP = logprob(cand) - logprob(curr)
  accept = logP > logU
  draw = accept ? cand : curr
  return (draw, accept)
end

end
