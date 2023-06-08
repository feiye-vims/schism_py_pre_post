from schism_py_pre_post.Grid.SourceSinkIn import source_sink, SourceSinkIn
from schism_py_pre_post.Grid.SMS import lonlat2cpp
from schism_py_pre_post.Timeseries.TimeHistory import TimeHistory
from schism_py_pre_post.Grid.Hgrid_extended import read_schism_hgrid_cached
import numpy as np
from scipy import spatial
import pickle


def nearest_neighbour(points_a, points_b):
    tree = spatial.cKDTree(points_b)
    return np.array(tree.query(points_a)[1]).reshape(-1,), np.array(tree.query(points_a)[0]).reshape(-1,)

def dist(points_group_A, points_group_B):
    points_A = np.squeeze(points_group_A.view(np.complex128))
    points_B = np.squeeze(points_group_B.view(np.complex128))
    return np.absolute(points_A-points_B)


def relocate_sources(
    old_ss_dir=None,
    feeder_info_file=None,
    hgrid_fname=None,
    outdir=None,
    max_search_radius = 2000.0,
    mandatory_sources_coor = np.empty((0, 2)),
    relocate_map = None,
    #--------------------------- end inputs -------------------------
):

    # Read original source/sink based on NWM
    old_source_sink = source_sink(source_dir=old_ss_dir)

    if relocate_map is None:
        # get old hgrid (without feeders)
        old_gd = read_schism_hgrid_cached(f'{old_ss_dir}/hgrid.gr3', overwrite_cache=False)

        # get old source coordinates
        old_gd.compute_ctr()

        # Note: source_eles starts from 1
        old_sources_coor = np.c_[old_gd.xctr[old_source_sink.source_eles-1], old_gd.yctr[old_source_sink.source_eles-1],
                                old_source_sink.vsource.get_time_average([])]
        np.savetxt(f'{old_ss_dir}/sources.xy', old_sources_coor)

        # Read feeder channel info
        with open(feeder_info_file, 'rb') as file:
            [feeder_l2g, feeder_points, feeder_heads, feeder_bases] = pickle.load(file)

        # get new hgrid (with feeders)
        new_gd = read_schism_hgrid_cached(hgrid_fname, overwrite_cache=False)
        new_gd.compute_ctr()

        # find matching source point at mandatory_sources_coor and feeders' base points
        new_sources_coor = np.r_[mandatory_sources_coor, feeder_bases[:, :2]]

        # get new source location (ele id)
        new_sources_eleids, _ = nearest_neighbour(np.r_[mandatory_sources_coor, feeder_heads[:, :2]], np.c_[new_gd.xctr, new_gd.yctr])

        # link each new source to the closest old source
        old_sources_x, old_sources_y = lonlat2cpp(old_sources_coor[:, 0], lat=old_sources_coor[:, 1])
        new_sources_x, new_sources_y = lonlat2cpp(new_sources_coor[:, 0], lat=new_sources_coor[:, 1])
        new2old_sources, relocation_distance = nearest_neighbour(
            np.c_[new_sources_x, new_sources_y],
            np.c_[old_sources_x, old_sources_y],
        )

        is_mandatory = np.zeros(len(new_sources_x, ), dtype=bool); is_mandatory[:len(mandatory_sources_coor)] = True
        nearby_relocation = np.logical_or(is_mandatory, relocation_distance < max_search_radius)

        new2old_sources[~nearby_relocation] = -1

        for old_source in np.unique(new2old_sources):
            if old_source == -1: continue  # skip far away sources

            ids = np.argwhere(new2old_sources == old_source)  # find all new sources mapped to the same old source
            if len(ids) == 0:
                print(f'old source {old_source} cannot be mapped to a new source')
            elif len(ids) == 1:  # exact match
                pass
            else:  # multiple new sources mapped to the same old source, pick the closest new source
                min_dist_id = np.argmin(relocation_distance[ids])
                new2old_sources[ids] = -1
                new2old_sources[ids[min_dist_id]] = old_source

        valid_relocation = new2old_sources >= 0
        valid_new_sources_eleids = new_sources_eleids[valid_relocation] + 1
        valid_new2old_sources = new2old_sources[valid_relocation]
        np.savetxt(f'{outdir}/relocate_map.txt', np.c_[valid_new_sources_eleids, valid_new2old_sources], fmt='%d %d')
    else:
        valid_new_sources_eleids = relocate_map[:, 0]; valid_new2old_sources = relocate_map[:, 1]
    
    # Assemble new source/sink files
    nsources = len(valid_new_sources_eleids)

    source_sink_in = SourceSinkIn(filename=None, number_of_groups=2, ele_groups=[valid_new_sources_eleids.tolist(),[]])
    vsource = TimeHistory(
        file_name=None,
        data_array=np.c_[old_source_sink.vsource.time, old_source_sink.vsource.data[:, valid_new2old_sources]],
        columns=['datetime'] + valid_new_sources_eleids.astype('str').tolist()
    )
    msource = TimeHistory(
        file_name=None,
        data_array=np.c_[
            old_source_sink.msource.time,
            -9999*np.ones([old_source_sink.msource.n_time, nsources]),
            np.zeros([old_source_sink.msource.n_time, nsources])
        ],
        columns=['datetime'] + valid_new_sources_eleids.astype('str').tolist() + valid_new_sources_eleids.astype('str').tolist()
    )

    source_sink_in.writer(f'{outdir}/source_sink.in')
    vsource.writer(f'{outdir}/vsource.th')
    msource.writer(f'{outdir}/msource.th')

if __name__ == "__main__":
    #--------------------------- inputs -------------------------
    old_ss_dir = '/sciclone/schism10/feiye/STOFS3D-v6/Inputs/v6_shadow_fcst/Relocate_SourceSink/original_source_sink/'
    feeder_info_file = f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/SMS_proj/feeder/feeder.pkl'
    hgrid_fname = '/sciclone/schism10/feiye/STOFS3D-v6/Inputs/I13v/Hgrid/hgrid.ll'
    outdir = '/sciclone/schism10/feiye/STOFS3D-v6/Inputs/v6_shadow_fcst/Relocate_SourceSink3/relocated_source_sink/'
    max_search_radius = 20000.0

    # Some major rivers may not have a feeder channel, or the feeder channel doesn't match the main river inputs.
    # In such cases, specify the actual river sources according to the original source locations
    # (e.g., by visualizing the map generated by see viz_source_.py).
    # These mandatory_sources_coor (along with the feeder channel bases) will be set as the new source locations
    mandatory_sources_coor = np.array([
        [-69.77256, 44.31494],  # Kennebec River, ME
        [-73.90869, 42.13509],  # Hudson River
        [-76.12903, 39.60946],  # Susquehanna River, VA
        [-74.94442, 40.34478],  # Delaware River
        [-78.425288, 34.508177],  # Cape Fear River, NC
        [-91.72306, 31.04462],  # Red River (upstream of Atchafalaya River)
        [-80.10808, 33.50005],  # Santee River
        [-79.81703, 33.59694],  # Black River
        [-79.57210, 33.71223],  # Black Mingo Creek
        [-79.49997, 33.84686],  # Lynches River
        [-79.48467, 33.93939],  # Pee Dee River
        [-79.33247, 33.98196],  # Little Pee Dee River
        [-77.917829, 34.749979],  # Northeast Cape Fear River
        [-87.9523, 30.8472]  # Mobile River
    ]).reshape(-1, 2)

    relocate_map = None

    relocate_sources(
        old_ss_dir=old_ss_dir,
        feeder_info_file=feeder_info_file,
        hgrid_fname=hgrid_fname,
        outdir=outdir,
        max_search_radius=max_search_radius,
        mandatory_sources_coor=mandatory_sources_coor,
        relocate_map=relocate_map
    )