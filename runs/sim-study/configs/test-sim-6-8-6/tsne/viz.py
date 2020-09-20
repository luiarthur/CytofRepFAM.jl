import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import sys
sys.path.append('../../../../PlotUtils')

def graph_tsne(tsne_df, clust, i, method, outpath, method_suffix=''):
    if clust is None:
        tsne_df[method] = tsne_df[method].astype(int)
    else:
        tsne_df[method] = clust.astype(int)
    mask_i = (tsne_df.sample_ind == i)
    df = tsne_df[mask_i]
    markersize = 15
    sns.pairplot(x_vars="tsne1", y_vars="tsne2", data=df,
                 hue=method,
                 plot_kws=dict(linewidth=0, s=markersize),
                 aspect=1, height=3)
    plt.savefig(outpath, bbox_inches="tight")
    plt.close()
    
def get_data_dict(tsne_df, y_header='Y', marker_header='M',
                  sample_ind_name='sample_ind',
                  true_labels_name='true_labels'):
    y_columns = filter(lambda x: x.startswith(y_header), tsne_df.columns)
    Y = tsne_df[y_columns].to_numpy()
    m_columns = filter(lambda x: x.startswith(marker_header), tsne_df.columns)
    M = tsne_df[m_columns].to_numpy() == 1
    Y[M] = np.nan
    sample_ind = tsne_df[sample_ind_name].to_numpy().astype(int)
    true_labels = tsne_df[true_labels_name].to_numpy().astype(int)
    return dict(Y=Y, M=M, sample_ind=sample_ind, true_labels=true_labels)

# TODO:
# make heatmaps
if __name__ == "__main__":
    path_to_csv = 'viz/csv'
    # methods = ['mclust', 'flowsom', 'rfam']
    methods = ['mclust', 'flowsom', 'true_labels']
    os.makedirs('viz/img', exist_ok=True)

    method = methods[0]

    for pmiss in [0.0, 0.2]:
        for zind in [1, 2, 3]:
            simname = f'pmiss{pmiss}-phi0-zind{zind}'

            tsne_path = f'{path_to_csv}/tsne-{simname}.csv'
            tsne_df = pd.read_csv(tsne_path)
            data_dict = get_data_dict(tsne_df)

            for method in methods:
                if method == 'true_labels':
                    clust = None
                else:
                    clust_path = f'{path_to_csv}/{method}-{simname}.csv'
                    print(clust_path)
                    clust = np.loadtxt(clust_path)
                for i in range(2):
                    outpath = f'viz/img/tsne-{method}{i + 1}-{simname}.pdf'
                    graph_tsne(tsne_df, clust, i + 1, method, outpath)

            # rfam
            method = "rfam"
            for phi in [0, 1, 10, 25, 100]:  # NOTE: mind this!
                clust_simname = f'pmiss{pmiss}-phi{phi}-zind{zind}'
                clust_path = f'{path_to_csv}/{method}-{clust_simname}.csv'
                print(clust_path)
                clust = np.loadtxt(clust_path)

                tsne_path = f'{path_to_csv}/tsne-{simname}.csv'
                tsne_df = pd.read_csv(tsne_path)
                for i in range(2):
                    outpath = f'viz/img/tsne-{method}{i + 1}-{clust_simname}.pdf'
                    graph_tsne(tsne_df, clust, i + 1, method, outpath,
                               method_suffix=f'phi={phi}')

