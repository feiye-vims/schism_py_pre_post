import numpy as np
from scipy.spatial import cKDTree


def inverse_distance_weighting(X, val, Xq, n_nei):
    my_ckdtree = cKDTree(X, leafsize=10)
    dist, idx = my_ckdtree.query(Xq, n_nei, eps=1e-6, p=2)
    dist[dist < 1e-9] = 1e-9  # avoid division by 0

    v = val[idx.ravel()].reshape(idx.shape)
    Vq = np.sum(v / dist, axis=1) / np.sum(1. / dist, axis=1)

    return Vq
