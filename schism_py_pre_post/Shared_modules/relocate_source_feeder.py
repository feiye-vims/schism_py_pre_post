from schism_py_pre_post.Grid.SourceSinkIn import source_sink, SourceSinkIn
from schism_py_pre_post.Grid.SMS import lonlat2cpp
from schism_py_pre_post.Timeseries.TimeHistory import TimeHistory
from schism_py_pre_post.Grid.Hgrid_ported import read_schism_hgrid_cached
import os
import numpy as np
from scipy import spatial
import pickle
import copy

from sqlalchemy import over


def nearest_neighbour(points_a, points_b):
    tree = spatial.cKDTree(points_b)
    return np.array(tree.query(points_a)[1]).reshape(-1,), np.array(tree.query(points_a)[0]).reshape(-1,)

def dist(points_group_A, points_group_B):
    points_A = np.squeeze(points_group_A.view(np.complex128))
    points_B = np.squeeze(points_group_B.view(np.complex128))
    return np.absolute(points_A-points_B)
    

if __name__ == "__main__":

    old_ss_dir = '/sciclone/schism10/feiye/STOFS3D-v6/Inputs/I23i/NWM/'

    # Read original source/sink based on NWM
    old_source_sink = source_sink(source_dir=old_ss_dir)

    # get old hgrid (without feeders)
    old_gd = read_schism_hgrid_cached(f'{old_ss_dir}/hgrid.gr3', overwrite_cache=False)

    # get old source coordinates
    old_gd.compute_ctr()
    old_sources_coor = np.c_[old_gd.xctr[old_source_sink.source_eles], old_gd.yctr[old_source_sink.source_eles], old_source_sink.vsource.get_time_average([])]
    np.savetxt(f'{old_ss_dir}/sources.xy', old_sources_coor)

    # Read feeder channel info
    with open(f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/SMS_proj/feeder/feeder.pkl', 'rb') as file:
        [feeder_l2g, feeder_points, feeder_heads, feeder_bases] = pickle.load(file)

    # get new hgrid (with feeders)
    new_gd = read_schism_hgrid_cached('/sciclone/schism10/feiye/STOFS3D-v6/Inputs/I23m/Hgrid/hgrid.ll', overwrite_cache=False)
    new_gd.compute_ctr()

    outdir = '/sciclone/schism10/feiye/STOFS3D-v6/Inputs/I23m/SourceSink_relocate/'

    # find matching source point at feeders' base points
    mandatory_sources_coor = np.array([
        [-69.77256, 44.31494],  # Kennebec River, ME
        [-73.90869, 42.13509],  # Hudson River
        [-74.94442, 40.34478],  # Delaware River
        [-78.425288, 34.508177],  # Cape Fear River, NC
        [-91.72306, 31.04462],  # Red River (upstream of Atchafalaya River)
        [-80.10808, 33.50005],  # Santee River
        [-79.81703, 33.59694],  # Black River
        [-79.57210, 33.71223],  # Black Mingo Creek
        [-79.49997, 33.84686],  # Lynches River
        [-79.48467, 33.93939],  # Pee Dee River
        [-79.33247, 33.98196],  # Little Pee Dee River
        [-77.917829, 34.749979]  # Northeast Cape Fear River
    ]).reshape(-1, 2)
    new_sources_coor = np.r_[mandatory_sources_coor, feeder_bases[:, :2]]

    # get new source location (ele id)
    new_sources_eleids, _ = nearest_neighbour(np.r_[mandatory_sources_coor, feeder_heads[:, :2]], np.c_[new_gd.xctr, new_gd.yctr])

    old_sources_x, old_sources_y = lonlat2cpp(old_sources_coor[:, 0], lat=old_sources_coor[:, 1])
    new_sources_x, new_sources_y = lonlat2cpp(new_sources_coor[:, 0], lat=new_sources_coor[:, 1])

    new2old_sources, relocation_distance = nearest_neighbour(
        np.c_[new_sources_x, new_sources_y],
        np.c_[old_sources_x, old_sources_y],
    )

    is_mandatory = np.zeros(len(new_sources_x, ), dtype=bool); is_mandatory[:len(mandatory_sources_coor)] = True
    nearby_relocation = np.logical_or(is_mandatory, relocation_distance < 2000.0)

    new2old_sources[~nearby_relocation] = -1

    for old_source in np.unique(new2old_sources):
        if old_source == -1: continue

        ids = np.argwhere(new2old_sources == old_source)
        if len(ids) == 0:
            print(f'old source {old_source} cannot be mapped to a new source')
        elif len(ids) == 1:
            pass
        else:
            min_dist_id = np.argmin(relocation_distance[ids])
            new2old_sources[ids] = -1
            new2old_sources[ids[min_dist_id]] = old_source

    # Assemble new source/sink files
    valid_relocation = new2old_sources >= 0
    nsources = sum(valid_relocation)
    eleids = new_sources_eleids[valid_relocation] + 1

    source_sink_in = SourceSinkIn(filename=None, number_of_groups=2, ele_groups=[eleids.tolist(),[]])
    vsource = TimeHistory(
        file_name=None,
        data_array=np.c_[old_source_sink.vsource.time, old_source_sink.vsource.data[:, new2old_sources[valid_relocation]]],
        columns=['datetime'] + eleids.astype('str').tolist()
    )
    msource = TimeHistory(
        file_name=None,
        data_array=np.c_[
            old_source_sink.msource.time,
            -9999*np.ones([old_source_sink.msource.n_time, nsources]),
            np.zeros([old_source_sink.msource.n_time, nsources])
        ],
        columns=['datetime'] + eleids.astype('str').tolist() + eleids.astype('str').tolist()
    )

    source_sink_in.writer(f'{outdir}/source_sink.in')
    vsource.writer(f'{outdir}/vsource.th')
    msource.writer(f'{outdir}/msource.th')

