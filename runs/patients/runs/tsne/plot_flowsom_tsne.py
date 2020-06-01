import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

path_to_flowsom = 'img/flowsom-clusterings.csv'
path_to_tsne = ['img/tsne-{}.txt'.format(i) for i in (1, 2, 3)]
flowsom = np.loadtxt(path_to_flowsom).astype(int)
fs = flowsom[:, 0]
sample_idx = flowsom[:, 1]

# TSNE for sampels combined
path_to_tsne_combined = 'img/tsne-combined.txt'
tsne = np.loadtxt(path_to_tsne_combined, delimiter=',')

# number of samples
num_sampels = len(path_to_tsne)

for i in range(num_sampels):
    print('Making figure {}'.format(i + 1))
    #
    mask_i = (sample_idx == i + 1)
    df = pd.DataFrame(tsne[mask_i], columns=["comp1", "comp2"])
    df["flowsom"] = fs[mask_i]
    # 
    markersize = 20 if df.shape[0] < 3000 else 10
    sns.pairplot(x_vars="comp1", y_vars="comp2", data=df,
                 hue="flowsom",
                 plot_kws=dict(linewidth=0, s=markersize),
                 aspect=1, height=5)
    plt.savefig("img/flowsom-tsne-combined-{}.pdf".format(i + 1),
                bbox_inches="tight")
    plt.close();

