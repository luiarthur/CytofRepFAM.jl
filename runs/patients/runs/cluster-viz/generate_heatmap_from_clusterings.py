import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sort_clusters
import sys

# Append path
sys.path.append("../../../PlotUtils")
import plot_yz

# Data dir
datadir = 'img'

# Path to data
path_to_data = ['img/y{}-prepped-with-missing.csv'.format(i)
                for i in (1, 2, 3)]

# Flowsom clusterings
fs_path = 'img/flowsom-clusterings.csv'
fs = np.loadtxt(fs_path).astype(int)

# Mclust clusterings
mclust_path = 'img/mclust-clusterings.csv'
mclust = np.loadtxt(mclust_path).astype(int)

# Read data
y = [pd.read_csv(path) for path in path_to_data]

# Number of cells in each sample
N = [yi.shape[0] for yi in y]

# Number of samples
I = len(y)

# Stack matrices
Y = np.vstack(y)

# sample indices
sample_idx = np.concatenate([np.full(N[i], i + 1) for i in range(I)])

def make_mclust_heatmap(i):
    print("plotting heatmap for sample {}".format(i + 1))
    idx_i = np.where(mclust[:, 1] == i + 1)[0]
    Yi = Y[idx_i, :]
    mclusti = mclust[idx_i, 0]
    wi = sort_clusters.compute_wi(mclusti)
    plot_yz.plot_y(Yi, wi, mclusti,
                   cm=plot_yz.blue2red.cm(9), vlim=(-4,4))
    plt.savefig('img/mclust-heatmap-{}.pdf'.format(i + 1), bbox="tight")
    plt.close()


def make_flowsom_heatmap(i):
    print("plotting heatmap for sample {}".format(i + 1))
    idx_i = np.where(fs[:, 1] == i + 1)[0]
    Yi = Y[idx_i, :]
    fsi = fs[idx_i, 0]
    wi = sort_clusters.compute_wi(fsi)
    plot_yz.plot_y(Yi, wi, fsi,
                   cm=plot_yz.blue2red.cm(9), vlim=(-4,4))
    plt.savefig('img/flowsom-heatmap-{}.pdf'.format(i + 1), bbox="tight")
    plt.close()


for i in range(I):
    make_mclust_heatmap(i)
    make_flowsom_heatmap(i)

