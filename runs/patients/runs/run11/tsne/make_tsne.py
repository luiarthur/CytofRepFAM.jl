import numpy as np
import pandas as pd
from sklearn.manifold import TSNE
from sklearn.preprocessing import scale
import os
import sys
from collections import Counter

sys.path.append('../../../../PlotUtils')
from Population import Population

def get_y_est(path_to_experiment, I):
    path = [f'{path_to_experiment}/img/txt/y{i}_mean.csv'
            for i in range(1, I + 1)]
    y = [np.loadtxt(p, delimiter=',') for p in path]
    return y

def compute_combined_tsne(y, seed=0):
    # Stack y for each sample into one matrix
    Y = np.vstack(y)

    # Compute tsne for big matrix (scaled)
    tsne = TSNE(verbose=2, random_state=seed).fit(scale(Y))

    # indices
    idx = [np.full(y[i].shape[0], i + 1) for i in range(len(y))]
    idx = np.concatenate(idx, axis=0)
    N = idx.shape[0]

    # First columns are Y, second and third to last are the embeddings, last
    # column are indices.
    output = np.concatenate((Y, tsne.embedding_, idx.reshape(-1, 1)),
                            axis=-1)

    # Create column header
    ncol = Y.shape[1]
    columns = [f'marker{j}' for j in np.arange(1, ncol + 1)]
    columns += ['tsne1', 'tsne2', 'sample_idx']

    return pd.DataFrame(output, columns=columns)

def relabel_lam(lams, Zs):
    assert len(lams) == len(Zs)
    I = len(lams)

    # Create a population dict.
    population = Population()

    # First, label most predominant subpopulations.
    for i in range(I):
        c = Counter(lams[i])
        for ci in c.most_common():
            k = ci[0] - 1
            population.label(Zs[i][:, k])

    # Now relabel after the most predominant subpopulations
    new_labels = [[population.label(Zs[i][:, lin - 1]) for lin in lams[i]]
                  for i in range(I)]

    return new_labels


# NOTE: Hardcoded...
def make_cluster_df(tsne_df, path_to_experiment):
    lam1 = np.loadtxt(f'{path_to_experiment}/img/txt/lam1_best.txt').astype(int)
    lam2 = np.loadtxt(f'{path_to_experiment}/img/txt/lam2_best.txt').astype(int)
    Z1 = np.loadtxt(f'{path_to_experiment}/img/txt/Z1_best.txt')
    Z2 = np.loadtxt(f'{path_to_experiment}/img/txt/Z2_best.txt')

    new_labels = relabel_lam(lams=[lam1, lam2], Zs=[Z1, Z2])
    new_labels = np.concatenate(new_labels)

    cluster_df = tsne_df.copy()
    cluster_df['cluster'] = new_labels

    return cluster_df


if __name__ == '__main__':
    nsamples = 2
    path_to_default_experiment = 'results/phi25-pi0.2'
    y = get_y_est(path_to_default_experiment, I=nsamples)
    tsne_df = compute_combined_tsne(y, seed=0)

    # Save results.
    os.makedirs('results/tsne', exist_ok=True)
    tsne_df.round(3).to_csv('results/tsne/tsne.csv')
