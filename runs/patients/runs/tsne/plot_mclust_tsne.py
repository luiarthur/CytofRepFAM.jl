import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os

os.makedirs('img/tsne', exist_ok=True)

path_to_mclust = ['img/mclust-{}.csv'.format(i) for i in (1, 2, 3)]
path_to_tsne = ['img/tsne-{}.txt'.format(i) for i in (1, 2, 3)]

for i in range(len(path_to_tsne)):
    print('Making figure {}'.format(i + 1))
    mclust = np.loadtxt(path_to_mclust[i]).astype(int)
    tsne = np.loadtxt(path_to_tsne[i], delimiter=',')
    df = pd.DataFrame(tsne, columns=["comp1", "comp2"])
    df["mclust"] = mclust
    
    markersize = 20 if df.shape[0] < 3000 else 10
    sns.pairplot(x_vars="comp1", y_vars="comp2", data=df, hue="mclust",
                 plot_kws=dict(linewidth=0, s=markersize),
                 aspect=1, height=5)
    plt.savefig("img/mclust-tsne-{}.pdf".format(i + 1), bbox_inches="tight")
    plt.close();

# TSNE for sampels combined
path_to_tsne_combined = 'img/tsne-combined.txt'
tsne = np.loadtxt(path_to_tsne_combined, delimiter=',')

n_lower = 0
n_upper = 0
for i in range(len(path_to_tsne)):
    print('Making figure {}'.format(i + 1))
    mclust = np.loadtxt(path_to_mclust[i]).astype(int)
    n_lower = n_upper
    n_upper = n_lower + mclust.shape[0]
    #
    df = pd.DataFrame(tsne[n_lower:n_upper], columns=["comp1", "comp2"])
    df["mclust"] = mclust
    # 
    markersize = 20 if df.shape[0] < 3000 else 10
    sns.pairplot(x_vars="comp1", y_vars="comp2", data=df, hue="mclust",
                 plot_kws=dict(linewidth=0, s=markersize),
                 aspect=1, height=5)
    plt.savefig("img/tsne/mclust-tsne-combined-{}.pdf".format(i + 1), bbox_inches="tight")
    plt.close();

