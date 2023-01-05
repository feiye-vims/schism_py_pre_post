from schism_py_pre_post.Grid.Hgrid_ported import read_schism_hgrid_cached
from schism_py_pre_post.Grid.SMS import get_all_points_from_shp
import pickle
from scipy import spatial
import numpy as np
from pathlib import Path
import os


def nearest_neighbour(points_a, points_b):
    tree = spatial.cKDTree(points_b)
    return np.array(tree.query(points_a)[1]).reshape(-1,), np.array(tree.query(points_a)[0]).reshape(-1,)

def dist(points_group_A, points_group_B):
    points_A = np.squeeze(points_group_A.view(np.complex128))
    points_B = np.squeeze(points_group_B.view(np.complex128))
    return np.absolute(points_A-points_B)

if __name__ == "__main__":
    # Read channel_centerline info
    channel_pts_xyz, _, _, _ = get_all_points_from_shp(fname='/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/SMS_proj/v14.42/auto_watershed_arcs+feeder_clipped_lonlat.shp')

    # get new hgrid (with feeders)
    gd = read_schism_hgrid_cached('/sciclone/schism10/feiye/STOFS3D-v6/Inputs/V6_mesh_from_JZ/Drag/drag.gr3', overwrite_cache=False)
    # pair channel points to grid points
    f2g, dist = nearest_neighbour(np.c_[channel_pts_xyz[:, 0], channel_pts_xyz[:, 1]], np.c_[gd.x, gd.y])

    # dredge channel points in hgrid
    gd.dp[f2g] = 0.0

    gd.save(os.path.splitext(gd.source_file)[0] + '.channel=0.gr3')

    pass
