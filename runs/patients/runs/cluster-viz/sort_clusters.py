import numpy as np
from collections import Counter

def order_by_cluster(clust):
    counter = Counter(clust)
    num_clusters = len(counter)
    c_dict = counter.most_common(num_clusters)
    cluster_labels_sorted = [k for (k, v) in c_dict]
    order = []
    for k in cluster_labels_sorted:
        order += np.where(clust == k)[0].tolist()
    return order

def compute_wi(lami):
    lami = np.array(lami)
    K = np.max(lami)
    wi = [np.mean(lami == k + 1) for k in range(K)]
    return np.array(wi)

