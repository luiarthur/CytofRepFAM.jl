import numpy as np
import pandas as pd
from sklearn.manifold import TSNE
from sklearn.preprocessing import scale

# Data dir
datadir = 'img'

# Path to data
path_to_data = ['img/y{}-prepped.csv'.format(i) for i in (1,2,3)]

# Read data
y = [pd.read_csv(path) for path in path_to_data]

# Fit tsne
tsne = [TSNE(verbose=2, random_state=0).fit(scale(yi)) for yi in y]

# Write results
for i in range(len(tsne)):
    np.savetxt('img/tsne-{}.txt'.format(i + 1),
               tsne[i].embedding_,
               fmt='%.5f', delimiter=',')

# Fit tsne
Y = np.vstack(y)
tsne_combined = TSNE(verbose=2, random_state=0).fit(scale(Y))

# Write results
np.savetxt('img/tsne-combined.txt',
           tsne_combined.embedding_,
           fmt='%.5f', delimiter=',')


