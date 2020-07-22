import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os

os.makedirs('img', exist_ok=True)

for phi in (0, 1, 100, 10000):  # (0, 1, 25, 50, 100, 1000)
    path_to_fam = 'img/fam-phi{}-clusterings.csv'.format(phi)
    fam_orig = np.loadtxt(path_to_fam).astype(int)
    fam = fam_orig[:, 0]
    sample_idx = fam_orig[:, 1]

    # TSNE for sampels combined
    path_to_tsne_combined = 'img/tsne-combined.txt'
    tsne = np.loadtxt(path_to_tsne_combined, delimiter=',')

    # number of samples
    num_sampels = np.unique(sample_idx).size

    for i in range(num_sampels):
        print('Making figure {}'.format(i + 1))
        #
        mask_i = (sample_idx == i + 1)
        df = pd.DataFrame(tsne[mask_i], columns=["comp1", "comp2"])
        clustername = "fam-phi{}".format(phi)
        df[clustername] = fam[mask_i]
        # 
        markersize = 20 if df.shape[0] < 3000 else 10
        sns.pairplot(x_vars="comp1", y_vars="comp2", data=df,
                     hue=clustername,
                     plot_kws=dict(linewidth=0, s=markersize),
                     aspect=1, height=5)
        plt.savefig("img/fam-phi{}-tsne-combined-{}.pdf".format(phi, i + 1),
                    bbox_inches="tight")
        plt.close();
