from schism_py_pre_post.Grid.Hgrid_ported import read_schism_hgrid_cached
from schism_py_pre_post.Grid.SMS import SMS_MAP
import pickle
from scipy import spatial
import numpy as np


def nearest_neighbour(points_a, points_b):
    tree = spatial.cKDTree(points_b)
    return np.array(tree.query(points_a)[1]).reshape(-1,), np.array(tree.query(points_a)[0]).reshape(-1,)

def dist(points_group_A, points_group_B):
    points_A = np.squeeze(points_group_A.view(np.complex128))
    points_B = np.squeeze(points_group_B.view(np.complex128))
    return np.absolute(points_A-points_B)

if __name__ == "__main__":
    # Read channel_arcs info
    center_lines = SMS_MAP(filename='/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/Outputs/CUDEM_merged_thalwegs_1e6_single_fix_simple_sms_cleaned_32cores/total_arcs.map')
    xyz, l2g = center_lines.get_xyz()

    # get old hgrid (without feeders)
    # no_feeder_gd = read_schism_hgrid_cached('/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/SMS_proj/no_feeder/hgrid.ll', overwrite_cache=False)

    # Read feeder channel info
    outdir = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/SMS_proj/feeder/'
    with open(f'{outdir}/feeder_arrays.pkl', 'rb') as file:
        [feeder_array_x, feeder_array_y] = pickle.load(file)

    feeder_centerline_xyz = np.zeros((0, 2), dtype=float)
    for x, y in zip(feeder_array_x, feeder_array_y):
        center_x = np.mean(x, axis=0)
        center_y = np.mean(y, axis=0)
        feeder_centerline_xyz = np.r_[feeder_centerline_xyz, np.c_[center_x, center_y]]
        pass

    # get new hgrid (with feeders)
    feeder_gd = read_schism_hgrid_cached('/sciclone/schism10/feiye/STOFS3D-v6/Inputs/V6_mesh_from_JZ/Test_levee_heights/hgrid.ll', overwrite_cache=False)

    ingrid = feeder_gd.inside_grid(xyz[:, :2]).astype(bool)
    centerline_xyz = np.r_[xyz[ingrid, :2], feeder_centerline_xyz]

    # pair channel points to grid points
    f2g, dist = nearest_neighbour(np.c_[centerline_xyz[:, 0], centerline_xyz[:, 1]], np.c_[feeder_gd.x, feeder_gd.y])

    # dredge channel points in hgrid
    feeder_gd.dp[f2g] += 1.0  # meters

    feeder_gd.save('/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/SMS_proj/feeder/hgrid_dredged.ll')

    pass

