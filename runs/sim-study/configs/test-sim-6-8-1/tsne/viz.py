import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# df = pd.DataFrame(np.random.randn(500, 2).cumsum(axis=0), columns=["x", "y"])
# sns.relplot(x="x", y="y", sort=False, kind="line", data=df);
# plt.show()

def graph_tsne(tsne_df, clust, i, method):
    NotImplemented
    
    
if __name__ == "__main__":
    path_to_csv = 'viz/csv'
    # methods = ['mclust', 'flowsom', 'rfam']
    methods = ['mclust', 'flowsom']

    tsne_df = '{}/pmiss0.6-phi0.0-zind1.csv'.format(path_to_csv)

