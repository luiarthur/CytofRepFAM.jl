import numpy as np
from sklearn.metrics import pairwise_distances
import matplotlib.pyplot as plt

def compare_xs(x1, x2, label=[None, None]):
    xmin = np.min([x1.min(), x2.min()])
    xmax = np.max([x1.max(), x2.max()])

    x1_counts = [(x1 == x).mean() for x in np.arange(xmin, xmax + 1)]
    plt.bar(np.arange(xmin, xmax+1) + 0.2, x1_counts, width=0.3, label=label[0])

    x2_counts = [(x2 == x).mean() for x in np.arange(xmin, xmax + 1)]
    plt.bar(np.arange(xmin, xmax+1) - 0.2, x2_counts, width=0.3, label=label[1])

    plt.xticks(np.arange(xmin, xmax+1))
    plt.legend()

def get_pairwise_dist(Z, W):
    _Z = Z[:, W > 0]
    J, K = _Z.shape
    z_pairwise_dist_mat = pairwise_distances(_Z.T, _Z.T, metric="l1")
    z_pairwise_dist = z_pairwise_dist_mat[np.triu_indices(K, k=1)]
    return z_pairwise_dist

def compute_Z_dist(Za, Zb, Wa, Wb, label=[None, None]):
    za_pairwise_dist = get_pairwise_dist(Za, Wa)
    zb_pairwise_dist = get_pairwise_dist(Zb, Wb)

    compare_xs(za_pairwise_dist, zb_pairwise_dist, label=label)
