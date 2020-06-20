import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append('../../../PlotUtils')
import plot_yz

data_dir = '../cluster-viz/img'
results_dir = 'tex/results/'

# Makes heatmaps of the small clusters

for i in [1, 2, 3]:
    for phi in [0, 1, 10, 25]:
        lami = np.loadtxt('{}/phi{}/img/txt/lam{}_best.txt'.format(results_dir, phi, i))
        wi = np.loadtxt('{}/phi{}/img/txt/W{}_best.txt'.format(results_dir, phi, i))
        K = wi.shape[0]
        yi = pd.read_csv('{}/y{}-prepped-with-missing.csv'.format(data_dir, i))
        df = pd.concat([yi, pd.DataFrame(dict(lami=lami))], axis=1)

        wi_df = pd.DataFrame(wi, index=np.arange(K) + 1)
        wi_df.columns=['wi']

        df = df.merge(wi_df, how='left', left_on='lami', right_index=True)

        # Largest cluster size to include
        w_upper = 0.05
        dfi = df[df.wi < w_upper]

        yi_small = dfi.values[:, :-2]
        lami_small = dfi.values[:, -2]

        wi[wi > w_upper] = 0

        plt.figure(figsize=(8, 8))
        plot_yz.plot_y(yi_small, wi, lami_small, vlim=(-4,4),
                       cm=plot_yz.blue2red.cm(7), interpolation='nearest')
        plt.savefig('out/phi{}-y{}.pdf'.format(phi, i))
        plt.close()
