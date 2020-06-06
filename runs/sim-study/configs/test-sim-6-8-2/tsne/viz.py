import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

def graph_tsne(tsne_df, clust, i, method, outpath, method_suffix=''):
    if clust is None:
        tsne_df[method] = tsne_df[method].astype(int)
    else:
        tsne_df[method] = clust.astype(int)
    mask_i = (tsne_df.sample_ind == i)
    df = tsne_df[mask_i]
    markersize = 20
    sns.pairplot(x_vars="tsne1", y_vars="tsne2", data=df,
                 hue=method,
                 plot_kws=dict(linewidth=0, s=markersize),
                 aspect=1, height=5)
    plt.savefig(outpath, bbox_inches="tight")
    plt.close()
    
    
if __name__ == "__main__":
    path_to_csv = 'viz/csv'
    # methods = ['mclust', 'flowsom', 'rfam']
    methods = ['mclust', 'flowsom', 'true_labels']
    os.makedirs('viz/img', exist_ok=True)

    method = methods[0]

    for pmiss in [0.0, 0.6]:
        for zind in [1, 2, 3]:
            simname = 'pmiss{}-phi0.0-zind{}'.format(pmiss, zind)
            for method in methods:
                if method == 'true_labels':
                    clust = None
                else:
                    clust_path = '{}/{}-{}.csv'.format(path_to_csv, method, simname)
                    print(clust_path)
                    clust = np.loadtxt(clust_path)
                tsne_path = '{}/tsne-{}.csv'.format(path_to_csv, simname)
                tsne_df = pd.read_csv(tsne_path)
                for i in range(2):
                    outpath = 'viz/img/tsne-{}{}-{}.pdf'.format(method, i + 1, simname)
                    graph_tsne(tsne_df, clust, i + 1, method, outpath)
            # rfam
            method = "rfam"
            for phi in [0.0, 1.0]:
                clust_simname = 'pmiss{}-phi{}-zind{}'.format(pmiss, phi, zind)
                clust_path = '{}/{}-{}.csv'.format(path_to_csv, method,
                                                   clust_simname)
                print(clust_path)
                clust = np.loadtxt(clust_path)

                tsne_path = '{}/tsne-{}.csv'.format(path_to_csv, simname)
                tsne_df = pd.read_csv(tsne_path)
                for i in range(2):
                    outpath = 'viz/img/tsne-{}{}-{}.pdf'.format(method,
                                                                i + 1,
                                                                clust_simname)
                    graph_tsne(tsne_df, clust, i + 1, method, outpath,
                               method_suffix='phi={}'.format(phi))

