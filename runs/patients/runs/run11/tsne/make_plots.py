# See: https://nikkimarinsek.com/blog/7-ways-to-label-a-cluster-plot-python
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = (5, 5)
from make_tsne import make_cluster_df


def plot_tsne(tsne_df, cluster_df, savepath, text_size=11):

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

def append_col(df, cluster):
    df = df.copy()
    df['cluster'] = cluster
    return df

if __name__ == '__main__':
    tsne_df = pd.read_csv('results/tsne/tsne.csv')
    text_size = 11

    # Plot FAM TSNE
    cluster_df = make_cluster_df(tsne_df, 'results/phi0')
    plot_tsne(tsne_df, cluster_df, 'results/tsne/tsne_rfam_phi0')

    # Plot repFAM TSNE
    cluster_df = make_cluster_df(tsne_df, 'results/phi25-pi0.2')
    plot_tsne(tsne_df, cluster_df, 'results/tsne/tsne_rfam_phi25')

    # Plot Mclust TSNE
    mclust = np.loadtxt('results/tsne/mclust.csv')[:, 0].astype(int)
    cluster_df = append_col(tsne_df, mclust)
    plot_tsne(tsne_df, cluster_df, 'results/tsne/tsne_mclust')

    # Plot FlowSOM TSNE
    flowsom = np.loadtxt('results/tsne/flowsom.csv')[:, 0].astype(int)
    cluster_df = append_col(tsne_df, mclust)
    plot_tsne(tsne_df, cluster_df, 'results/tsne/tsne_flowsom')