# NOTE: This script prints the cooccurence of missing expression levels
# for the CB data.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

path_to_cb_data = '../../runs/data/cb_transformed.csv'
cb = pd.read_csv(path_to_cb_data)
ms = [(g[1].drop(['sample_id'], axis=1).isna() * 1.0).to_numpy()
      for g in cb.groupby('sample_id')]

i = 0
mi_counts = ms[i].T.dot(ms[i])
plt.imshow(mi_counts / ms[i].shape[0])
plt.colorbar()
plt.show()
