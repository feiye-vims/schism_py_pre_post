from pylib import schism_grid  # from ZG's pylib: pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple pylibs4schism==0.1.10
from scipy import spatial
import numpy as np


def nearest_neighbour(points_a, points_b):
    tree = spatial.cKDTree(points_b)
    return tree.query(points_a)[1]


if __name__ == "__main__":
    xyz = np.loadtxt('/sciclone/schism10/feiye/STOFS3D-v4/AWIPS/AWIPS_mask_from_Yuji/conus.xyz')
    xyz[xyz[:, 0]>180, 0] -= 360
    grid = schism_grid('/sciclone/schism10/feiye/STOFS3D-v4/Inputs/v11.7/hgrid.ll.pkl')

    ie, ip, acor = grid.compute_acor(np.c_[xyz[:, 0], xyz[:, 1]])
    ip[ie==-1, :] *= -1

    fmt = '%12d '*5 + '%20.10f '*3
    np.savetxt('/sciclone/schism10/feiye/STOFS3D-v4/AWIPS/AWIPS_mask_from_Yuji/mask0.txt', np.c_[np.array(range(len(ie)))+1, ie, ip, acor], fmt=fmt)
    fmt = '%12d '*4 + '%20.10f '*3
    np.savetxt('/sciclone/schism10/feiye/STOFS3D-v4/AWIPS/AWIPS_mask_from_Yuji/mask.txt', np.c_[np.array(range(len(ie)))+1, ip, acor], fmt=fmt)
    
    pass