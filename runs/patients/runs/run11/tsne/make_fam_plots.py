# See: https://nikkimarinsek.com/blog/7-ways-to-label-a-cluster-plot-python
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = (5, 5)
from make_tsne import make_cluster_df

def plot_tsne(tsne_df, path_to_clustering, savepath, text_size=11):
    cluster_df = make_cluster_df(tsne_df, path_to_clustering)

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

        plt.savefig(f'{savepath}_sample{int(i)}.pdf')
        plt.close()

if __name__ == '__main__':
    tsne_df = pd.read_csv('results/tsne/tsne.csv')
    text_size = 11

    # Plot FAM TSNE
    plot_tsne(tsne_df, 'results/phi0',
              'results/tsne/tsne_rfam_phi0')

    # Plot repFAM TSNE
    plot_tsne(tsne_df, 'results/phi25-pi0.2',
              'results/tsne/tsne_rfam_phi25')

    # Plot Mclust TSNE
    # Plot FlowSOM TSNE
