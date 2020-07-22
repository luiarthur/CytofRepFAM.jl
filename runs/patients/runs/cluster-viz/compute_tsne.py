import numpy as np
import pandas as pd
from sklearn.manifold import TSNE
from sklearn.preprocessing import scale
import os

# Data dir
datadir = 'img'
# datadir = '{}/cytof/results/repfam/patients-data/run11/phi0/img/txt'.format(os.environ['SCRATCH_DIR'])

# Path to data
# path_to_data = ['img/y{}-prepped.csv'.format(i) for i in (1,2,3)]
path_to_data = ['{}/y{}-prepped.csv'.format(datadir, i) for i in (1,2)]
# path_to_data = ['{}/y{}_mean.csv'.format(datadir, i) for i in (1,2)]

# Read data
# y = [pd.read_csv(path, header=None) for path in path_to_data]
y = [pd.read_csv(path) for path in path_to_data]

# for i in range(len(y)):
#     y[i].round(3).to_csv('{}/y{}-prepped.csv'.format('img', i + 1))

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

# indices
idx = [np.full(y[i].shape[0], i) for i in range(len(y))]
idx = np.concatenate(idx, axis=0)
N = idx.shape[0]

# First columns are Y, second and third to last are the embeddings, last column
# are indices.
output = np.concatenate((Y, tsne_combined.embedding_, idx.reshape(N, 1)),
                        axis=-1)

# Write results
np.savetxt('img/tsne-combined.txt',
           tsne_combined.embedding_,
           fmt='%.5f', delimiter=',')

