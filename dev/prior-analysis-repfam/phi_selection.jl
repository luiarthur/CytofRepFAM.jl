import Pkg; Pkg.activate("../../")

import CytofRepFAM.Model.ImportanceSampling: exp_f_rfam, min_pairwise_dist_ge_u, rep_fn

phis = (0, 5, 10, 15, 20, 25, 30)
us = 0:5
exp_f = [[let
            println("phi: $phi, u: $u")
            exp_f_rfam(N=10000, J=15, K=25, 
                       repulsive_fn=Z->rep_fn(Z, phi),
                       f=Z->min_pairwise_dist_ge_u(Z, u))
          end for u in us] for phi in phis]

