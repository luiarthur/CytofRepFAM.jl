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

import os
import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('../../../../../PlotUtils')
import blue2red
import plot_yz
from Population import Population
import zinfo
import tqdm
from z_dist import compute_Z_dist


def genpaths(results_dir):
    return [f'{results_dir}/{p}' for p in os.listdir(results_dir)]


if __name__ == '__main__':
    results_dir = '../results'
    paths = genpaths(results_dir)
    paths = [path for path in paths if 'phi' in path]
    markernames=["2B4", "3DL1", "CD158B", "CD8", "CD94",
                 "CKIT", "DNAM1", "EOMES", "NKG2A", "NKG2D" ,
                 "NKP30", "SIGLEC7", "SYK", "TBET", "ZAP70"]

    # Compute Z dist
    os.makedirs(f'{results_dir}/misc', exist_ok=True)
    for i in (1, 2):
        patha = f'{results_dir}/phi25-pi0.2/img/txt'
        zai = np.loadtxt(f'{patha}/Z{i}_best.txt')
        wai = np.loadtxt(f'{patha}/W{i}_best.txt')

        pathb = f'{results_dir}/phi0/img/txt'
        zbi = np.loadtxt(f'{pathb}/Z{i}_best.txt')
        wbi = np.loadtxt(f'{pathb}/W{i}_best.txt')

        plt.figure(figsize=(5, 5))
        compute_Z_dist(zai, zbi, wai, wbi, label=['rep-FAM', 'ind-FAM'])
        plt.xlabel("pairwise distance ")
        plt.ylabel("probability")
        plt.savefig(f'{results_dir}/misc/Z{i}_dist_compare.pdf',
                    bbox_inches="tight")
        plt.close()

    # Make other graphs.
    for path in tqdm.tqdm(paths):
        population = Population()

        R_path = f'{path}/img/txt/Rcounts.txt'
        R = np.loadtxt(R_path)
        plt.figure(figsize=(5,5))
        zinfo.plot_num_selected_features(R, ymax=0.8, 
                                         xlabel='# of selected subpopulations')
        plt.savefig(f'{path}/img/Rcounts.pdf', bbox_inches="tight")
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
            Zi = np.loadtxt(zi_path)

            # Plot y centroids
            plt.figure(figsize=(6,6))
            plot_yz.plot_y_centroids(yi, lami, wi, vlim=(-3, 3), cm=blue2red.cm(6),
                                     population=population, Zi=Zi, 
                                     fs_xlabel=16, fs_ylabel=16,
                                     fs_xticks=14, fs_yticks=16,
                                     rotation=45, markernames=markernames,
                                     ha="right")
            outpath = f'{path}/img/y{i}_centroid.pdf'
            plt.savefig(outpath, bbox_inches="tight")
            plt.close()
