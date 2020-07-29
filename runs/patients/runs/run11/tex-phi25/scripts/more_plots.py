# We may make a figure similar to figures 4 and 5 but with means of y_inj in
# the cluster, instead of z_jk.  e.g., using the heatmaps in figure 6, we can
# compute the means of y for each (feature * marker), ybar_jk, marker j and
# feature k.  We then make a plot of ybar_jk similar to the plot in figure 4.
# Instead of 0 or 1, we have actual expression levels at some aggregated level.
# We can keep ordering of features according to their abundances, but we may
# let the R heatmap function re-arrange the markers to show the difference
# between the features (any order of the markers does not have any meaning).
#
# This figure can be away from having 0 or 1. Also, it can show all features
# well.  The current heatmaps of y do not show well features with low
# abundance. 

import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('../../../../../PlotUtils')
import blue2red
import plot_yz
from Population import Population

if __name__ == '__main__':
    results_path = '../results'

    for phi in (25, ):
        for pi in (0.0, 0.1, 0.2, 0.3):
            population = Population()
            for i in (1, 2):
                # Read data
                yi_path = '{}/phi{}-pi{}/img/txt/y{}_mean.csv'.format(results_path, phi, pi, i)
                yi = np.loadtxt(yi_path, delimiter=',')
                lami_path = '{}/phi{}-pi{}/img/txt/lam{}_best.txt'.format(results_path, phi, pi, i)
                lami = np.loadtxt(lami_path, dtype=int)
                wi_path = '{}/phi{}-pi{}/img/txt/W{}_best.txt'.format(results_path, phi, pi, i)
                wi = np.loadtxt(wi_path)
                zi_path = '{}/phi{}-pi{}/img/txt/Z{}_best.txt'.format(results_path, phi, pi, i)
                Zi = np.loadtxt(zi_path)

                # Plot
                plt.figure(figsize=(6,6))
                plot_yz.plot_y_centroids(yi, lami, wi, vlim=(-3, 3), cm=blue2red.cm(6),
                                         population=population, Zi=Zi, 
                                         fs_xlabel=16, fs_ylabel=16,
                                         fs_xticks=14, fs_yticks=16)
                outpath = '{}/phi{}-pi{}/img/y{}_centroid.pdf'.format(results_path, phi, pi, i)
                plt.savefig(outpath, bbox_inches="tight")
                plt.close()
