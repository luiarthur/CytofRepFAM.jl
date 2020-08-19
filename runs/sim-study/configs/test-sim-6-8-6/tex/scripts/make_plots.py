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
import os
import sys
import re

sys.path.append('../../../../../PlotUtils')
import blue2red
import plot_yz
from Population import Population
import zinfo
import tqdm

def genpaths(results_dir):
    return [f'{results_dir}/{p}' for p in os.listdir(results_dir)]

def get_zind(path):
    return re.findall(r'(?<=zind)\d+', path)[0]

if __name__ == '__main__':
    results_dir = '../results'
    true_Z_path = '../../'
    paths = [p for p in genpaths(results_dir) if 'pmiss' in p]

    for path in tqdm.tqdm(paths):
        # Create population indices.
        population = Population()
        zind = get_zind(path)
        trueZ = np.loadtxt(f'{true_Z_path}/Z{zind}.txt')
        for k in range(trueZ.shape[1]):
            _ = population.label(trueZ[:, k])

        # Rcounts
        R_path = f'{path}/img/txt/Rcounts.txt'
        R = np.loadtxt(R_path)
        plt.figure(figsize=(4,4))
        zinfo.plot_num_selected_features(R, refs=[6, 5], ymax=1.05)
        plt.savefig(f'{path}/img/Rcounts.pdf', bbox_inches='tight')
        plt.close()

        for i in (1, 2):
            # Read data
            yi_path = f'{path}/img/txt/y{i}_mean.csv'
            yi = np.loadtxt(yi_path, delimiter=',')
            lami_path = f'{path}/img/txt/lam{i}_best.txt'
            lami = np.loadtxt(lami_path, dtype=int)
            wi_path = f'{path}/img/txt/W{i}_best.txt'
            wi = np.loadtxt(wi_path)
            zi_path = f'{path}/img/txt/Z{i}_best.txt'
            zi = np.loadtxt(zi_path)

            # Plot
            plt.figure(figsize=(6,6))
            plot_yz.plot_y_centroids(yi, lami, wi, vlim=(-3, 3), cm=blue2red.cm(6),
                                     population=population, Zi=zi,
                                     fs_xlabel=16, fs_ylabel=16,
                                     fs_xticks=16, fs_yticks=16)
            outpath = f'{path}/img/y{i}_centroid.pdf'
            plt.savefig(outpath, bbox_inches="tight")
            plt.close()

            # Plot Z estimate.
            plt.figure(figsize=(6,6))
            plot_yz.plot_Z(Z_mean=zi, wi_mean=wi, lami_est=lami, w_thresh=0.0,
                           population=population, add_colorbar=False,
                           ylab='features (abundance)',
                           fs_lab=15, fs_celltypes=15, fs_markers=15,
                           fs_cbar=15)
            outpath = f'{path}/img/Z{i}.pdf'
            plt.savefig(outpath, bbox_inches="tight")
            plt.close()
