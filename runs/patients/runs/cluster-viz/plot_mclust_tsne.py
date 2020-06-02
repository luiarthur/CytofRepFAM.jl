import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os

os.makedirs('img', exist_ok=True)

path_to_mclust = 'img/mclust-clusterings.csv'
mclust = np.loadtxt(path_to_mclust).astype(int)
mc = mclust[:, 0]
sample_idx = mclust[:, 1]

# TSNE for samples combined
path_to_tsne_combined = 'img/tsne-combined.txt'
tsne = np.loadtxt(path_to_tsne_combined, delimiter=',')

# number of samples
num_samples = np.unique(sample_idx).size

for i in range(num_samples):
    print('Making figure {}'.format(i + 1))
    #
    mask_i = (sample_idx == i + 1)
    df = pd.DataFrame(tsne[mask_i], columns=["comp1", "comp2"])
    df["mclust"] = mc[mask_i]
    # 
    markersize = 20 if df.shape[0] < 3000 else 10
    sns.pairplot(x_vars="comp1", y_vars="comp2", data=df,
                 hue="mclust",
                 plot_kws=dict(linewidth=0, s=markersize),
                 aspect=1, height=5)
    plt.savefig("img/mclust-tsne-combined-{}.pdf".format(i + 1),
                bbox_inches="tight")
    plt.close();

