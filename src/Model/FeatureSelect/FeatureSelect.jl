include("DataFS.jl")
include("ConstantsFS.jl")
include("StateFS.jl")
include("TunersFS.jl")

include("update_p_feature_select.jl")
include("update_omega_feature_select.jl")
include("update_W_feature_select.jl")
include("update_r_feature_select.jl")
include("update_lam_feature_select.jl")

include("compute_marg_loglike.jl")
include("update_feature_select.jl")
include("fit_feature_select.jl")

include("swap_chains.jl")
include("fit_feature_select_pt.jl")

include("sample_minibatch.jl")
include("update_via_trained_prior.jl")
include("fit_feature_select_train_prior.jl")
include("fit_feature_select_imcmc_pt.jl")
