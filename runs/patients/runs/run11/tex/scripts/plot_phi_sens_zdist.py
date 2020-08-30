import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams["font.size"] = 15
mpl.rcParams["xtick.labelsize"] = 15
mpl.rcParams["ytick.labelsize"] = 15
mpl.rcParams["figure.figsize"] = (6, 6)

import seaborn as sns
import z_dist
from collections import Counter


def get_pairwise_dist(i, phi, results_dir):
    path = f'{results_dir}/phi{phi}/img/txt'
    zi = np.loadtxt(f'{path}/Z{i}_best.txt')
    wi = np.loadtxt(f'{path}/W{i}_best.txt')
    zdi = z_dist.get_pairwise_dist(zi, wi).astype(int)
    return pd.DataFrame(dict(distance=zdi, phi=phi, i=i))


def plot_phi_sens_zdist(i, phis, results_dir, I=2):
    dfs = [get_pairwise_dist(i, phi, results_dir) for phi in phis]
    df = pd.concat(dfs)
    df = (df.groupby(['phi'])['distance']
            .value_counts(normalize=True)
            .rename('probability')
            .reset_index()
            .sort_values('distance'))
    sns.barplot(x="distance", y="probability", hue="phi", data=df)
    plt.legend(loc="upper right", title="$\phi$")
    return df


if __name__ == '__main__':
    results_dir = '../results'
    phis = (1, 10, 100)
    I = 2
    for i in range(1, I + 1):
        df = plot_phi_sens_zdist(i, phis, results_dir)
        plt.xticks(rotation=90)
        plt.savefig(f'{results_dir}/misc/z{i}_dist_phi_sens.pdf',
                    bbox_inches="tight")
        plt.close()
