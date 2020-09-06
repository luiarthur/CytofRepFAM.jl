# See: https://nikkimarinsek.com/blog/7-ways-to-label-a-cluster-plot-python
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = (5, 5)
from make_tsne import make_cluster_df

sys.path.append('../../../../PlotUtils')
import plot_yz

def plot_tsne(tsne_df, cluster_df, savepath, text_size=11, plot_y=False, m=None):
    n_clus = len(cluster_df.cluster.unique())
    cp = sns.color_palette('Paired', n_colors=n_clus)

    for i in cluster_df.sample_idx.unique():
        # data = cluster_df[cluster_df.sample_idx == i].sample(1000)
        data = cluster_df[cluster_df.sample_idx == i]
        for _, label in enumerate(data.cluster.unique()):
            # loop through data points and plot each point 
            d = data[data.cluster == label]
            # add the data point as text
            plt.scatter(d.tsne1, d.tsne2, color=cp[label - 1], s=6)
            centroid = data.loc[data['cluster']==label,
                                ['tsne1','tsne2']].median()
            t = plt.annotate(int(label), centroid,
                             horizontalalignment='center',
                             verticalalignment='center',
                             size=text_size, weight='bold',
                             color='black', 
                             backgroundcolor=cp[label - 1])
            t.set_bbox(dict(facecolor=cp[label - 1], alpha=0.6, linewidth=0))

        plt.savefig(f'{savepath}_sample{int(i)}.pdf', bbox_inches="tight")
        plt.close()

        if plot_y:
            K = 25  # Hardcoded.
            colnames = [c for c in data.columns if 'marker' in c]
            yi = data[colnames].to_numpy() + 0

            if m is not None:
                mi = m[cluster_df.sample_idx == i]
                yi[mi == 1] = np.nan

            markernames=["2B4", "3DL1", "CD158B", "CD8", "CD94",
                         "CKIT", "DNAM1", "EOMES", "NKG2A", "NKG2D" ,
                         "NKP30", "SIGLEC7", "SYK", "TBET", "ZAP70"]
            lami_est = data.cluster

            plot_yz.plot_y(yi, K, lami_est, vlim=(-4, 4),
                           markernames=markernames, ha='right',
                           rotation=45, interpolation='nearest')
            plt.savefig(f'{savepath}_y{int(i)}.pdf', bbox_inches="tight")
            plt.close()



def append_col(df, cluster):
    df = df.copy()
    df['cluster'] = cluster
    return df

if __name__ == '__main__':
    tsne_df = pd.read_csv('results/tsne/tsne.csv')
    text_size = 11
    m = np.loadtxt('results/misc/m.txt', delimiter=',')

    # Plot FAM TSNE
    cluster_df = make_cluster_df(tsne_df, 'results/phi0')
    plot_tsne(tsne_df, cluster_df, 'results/tsne/tsne_rfam_phi0')

    # Plot repFAM TSNE
    cluster_df = make_cluster_df(tsne_df, 'results/phi25-pi0.2')
    plot_tsne(tsne_df, cluster_df, 'results/tsne/tsne_rfam_phi25')

    # Plot Mclust TSNE
    mclust = np.loadtxt('results/tsne/mclust.csv')[:, 0].astype(int)
    cluster_df = append_col(tsne_df, mclust)
    plot_tsne(tsne_df, cluster_df, 'results/tsne/tsne_mclust', plot_y=True, m=m)

    # Plot FlowSOM TSNE
    flowsom = np.loadtxt('results/tsne/flowsom.csv')[:, 0].astype(int)
    cluster_df = append_col(tsne_df, flowsom)
    plot_tsne(tsne_df, cluster_df, 'results/tsne/tsne_flowsom', plot_y=True, m=m)
