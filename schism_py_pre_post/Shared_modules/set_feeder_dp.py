from schism_py_pre_post.Grid.SourceSinkIn import source_sink, SourceSinkIn
from schism_py_pre_post.Grid.SMS import lonlat2cpp
from schism_py_pre_post.Timeseries.TimeHistory import TimeHistory
from schism_py_pre_post.Grid.Hgrid_ported import read_schism_hgrid_cached
import os
import numpy as np
from scipy import spatial
import pickle
import copy


def nearest_neighbour(points_a, points_b):
    tree = spatial.cKDTree(points_b)
    return np.array(tree.query(points_a)[1]).reshape(-1,), np.array(tree.query(points_a)[0]).reshape(-1,)

def dist(points_group_A, points_group_B):
    points_A = np.squeeze(points_group_A.view(np.complex128))
    points_B = np.squeeze(points_group_B.view(np.complex128))
    return np.absolute(points_A-points_B)

if __name__ == "__main__":
    # Read feeder channel info
    input_dir = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/SMS_proj/feeder/'
    with open(f'{input_dir}/feeder.pkl', 'rb') as file:
        [feeder_l2g, feeder_points, feeder_heads, feeder_bases] = pickle.load(file)

    # get old hgrid (without feeders)
    old_gd = read_schism_hgrid_cached('/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/SMS_proj/no_feeder/hgrid.ll', overwrite_cache=False)

    # get new hgrid (with feeders)
    gd = read_schism_hgrid_cached('/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/SMS_proj/feeder/hgrid.ll', overwrite_cache=False)

    # pair feeder channel points to grid points
    f2g, _ = nearest_neighbour(feeder_points, np.c_[gd.x, gd.y])

    # find outside grid points
    feeder_in_grid = old_gd.inside_grid(feeder_points).astype(bool)

    for i, id in enumerate(feeder_l2g):
        gd_points_in_feeder = f2g[id]
        gd_points_in_feeder_in_grid = feeder_in_grid[id]
        gd_points_in_external_feeder = gd_points_in_feeder[~gd_points_in_feeder_in_grid]

        gd.dp[gd_points_in_external_feeder] = np.min(gd.dp[gd_points_in_feeder])
    
    gd.save(f'{input_dir}/hgrid_feeder_dp_set.ll')
    pass


