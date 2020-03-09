include("imports.jl")

path_to_output = ("/scratchdata/alui2/cytof/results/repfam/test-sim-6-7-2/" *
                  "maxtemp10-ntempts20-degree2-N500/output.bson")

output = BSON.load(path_to_output)

# 1. For a fixed state (theta), compute_marg_loglike for different temperatures.
tempers = output[:tempers]
# tempers = CytofRepFAM.MCMC.WSPT.gentempers(2, 20, degree=1)
lls = [CytofRepFAM.Model.compute_marg_loglike(output[:lastState].theta,
                                              output[:c].constants,
                                              output[:d].data,
                                              tau) for tau in tempers];
plt.plot(tempers, lls, marker="o")
plt.savefig("tmp.pdf", bbox_inches="tight")
plt.close()

# plt.plot(log.(tempers), lls - [ll[end] for ll in output[:lls]], marker="o")
# plt.savefig("tmp.pdf", bbox_inches="tight")
# plt.close()
