# See: https://nikkimarinsek.com/blog/7-ways-to-label-a-cluster-plot-python
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = (5, 5)
from make_tsne import make_cluster_df

path_to_experiments = ['results/phi0', 'results/phi25-pi0.2']
tsne_df = pd.read_csv('results/tsne/tsne.csv')
text_size = 11

# Plot FAM TSNE
for path in path_to_experiments:
    cluster_df = make_cluster_df(tsne_df, path)

    # for i in cluster_df.sample_idx.unique():
    #     sns.pairplot(x_vars="tsne1", y_vars="tsne2",
    #                  data=cluster_df[cluster_df.sample_idx == i],
    #                  hue='cluster',
    #                  plot_kws=dict(linewidth=0), aspect=1, height=3)
    #     plt.tight_layout()
    #     plt.show()

    # loop through labels and plot each cluster

    # TODO: 
    # - Try different color pallete
    # - Try sampling the data to visualize easier
    n_clus = len(cluster_df.cluster.unique())
    cp = sns.color_palette('muted', n_colors=n_clus)

    for i in cluster_df.sample_idx.unique():
        data = cluster_df[cluster_df.sample_idx == i].sample(1000)
        for j, label in enumerate(data.cluster.unique()):
            # loop through data points and plot each point 
            d = data[data.cluster == label]
            # add the data point as text
            plt.scatter(d.tsne1, d.tsne2, color=cp[j], s=10) 
            centroid = data.loc[data['cluster']==label,
                                ['tsne1','tsne2']].mean()
            t = plt.annotate(int(label), centroid,
                             horizontalalignment='center',
                             verticalalignment='center',
                             size=text_size, weight='bold',
                             color='white', 
                             backgroundcolor=cp[j])
            t.set_bbox(dict(facecolor=cp[j], alpha=0.5, linewidth=0))

        # plt.xlim(cluster_df.tsne1.min(), cluster_df.tsne1.max())
        # plt.ylim(cluster_df.tsne2.min(), cluster_df.tsne2.max())
        plt.show()
