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
    mandatory_sources_coor = np.empty((0, 4)),
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
        np.savetxt(f'{old_ss_dir}/original_sources.xy', old_sources_coor)

        # Read feeder channel info
        with open(feeder_info_file, 'rb') as file:
            [feeder_l2g, feeder_points, feeder_heads, feeder_bases] = pickle.load(file)
        np.savetxt(f'{outdir}/feeder_heads.xy', feeder_heads)
        np.savetxt(f'{outdir}/feeder_bases.xy', feeder_bases)

        # get new hgrid (with feeders)
        new_gd = read_schism_hgrid_cached(hgrid_fname, overwrite_cache=False)
        new_gd.compute_ctr()

        # process nan values in mandatory_sources_coor
        for i, row in enumerate(mandatory_sources_coor):
            if np.isnan(row[2]):
                mandatory_sources_coor[i, 2] = mandatory_sources_coor[i, 0]
            if np.isnan(row[3]):
                mandatory_sources_coor[i, 3] = mandatory_sources_coor[i, 1]

        # find matching source point at mandatory_sources_coor and feeders' base points
        # these are the desired new source locations
        new_sources_coor = np.r_[mandatory_sources_coor[:, :2], feeder_heads[:, :2]]
        # these are the locations used to search for the closest old source
        # , in other words, to link the new source to the old source
        new_sources_search_point_coor = np.r_[mandatory_sources_coor[:, 2:4], feeder_bases[:, :2]]

        # get new source location (ele id, index starting from 0)
        new_sources_eleids, _ = nearest_neighbour(new_sources_coor, np.c_[new_gd.xctr, new_gd.yctr])

        # link each new source to the closest old source
        old_sources_x, old_sources_y = lonlat2cpp(old_sources_coor[:, 0], lat=old_sources_coor[:, 1])
        new_sources_x, new_sources_y = lonlat2cpp(new_sources_search_point_coor[:, 0], lat=new_sources_search_point_coor[:, 1])
        new2old_sources, relocation_distance = nearest_neighbour(
            np.c_[new_sources_x, new_sources_y],
            np.c_[old_sources_x, old_sources_y],
        )

        is_mandatory = np.zeros(len(new_sources_x, ), dtype=bool); is_mandatory[:len(mandatory_sources_coor)] = True
        # mandatory sources are always relocated despite the distance
        valid_relocation = np.logical_or(is_mandatory, relocation_distance < max_search_radius)
        # some new sources are not possible to link to old sources and marked as -1
        new2old_sources[~valid_relocation] = -1

        for old_source in np.unique(new2old_sources):
            if old_source == -1: continue  # skip invalid new sources

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
    # Note: source_eles starts from 1
    new_sources_coor = np.c_[
        new_gd.xctr[source_sink_in.ip_group[0]-1],
        new_gd.yctr[source_sink_in.ip_group[0]-1],
        vsource.get_time_average([])
    ]
    np.savetxt(f'{outdir}/relocated_sources.xyz', new_sources_coor)

    source_sink_in.writer(f'{outdir}/source_sink.in')
    vsource.writer(f'{outdir}/vsource.th')
    msource.writer(f'{outdir}/msource.th')

if __name__ == "__main__":
    #--------------------------- inputs -------------------------
    old_ss_dir = '/sciclone/schism10/feiye/STOFS3D-v6/Inputs/v6_shadow_fcst/Relocate_SourceSink/original_source_sink/'
    feeder_info_file = f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/SMS_proj/feeder/feeder.pkl'
    hgrid_fname = '/sciclone/schism10/feiye/STOFS3D-v6/Inputs/I13v/Hgrid/hgrid.ll'
    outdir = '/sciclone/schism10/feiye/STOFS3D-v6/Inputs/v6_shadow_fcst/Relocate_SourceSink3/relocated_source_sink/'
    max_search_radius = 2000.0

    # Some major rivers may not have a feeder channel, or the feeder channel doesn't match the main river inputs.
    # In such cases, manually relocate sources with the mandatory_sources_coor array.
    # (e.g., by visualizing the map generated by see viz_source.py).
    # The first two columns are longitude and latitude of the mandatory (relocated) source locations.
    # The third and fourth columns are longitude and latitude of the corresponding original source locations.
    # If the target original source locations are the same as the mandatory source locations,
    # the third and fourth columns can be set as np.nan
    mandatory_sources_coor = np.array([
        [-69.77256, 44.31494, np.nan, np.nan],  # Kennebec River, ME
        [-73.90869, 42.13509, np.nan, np.nan],  # Hudson River, NY
        [-76.12903, 39.60946, np.nan, np.nan],  # Susquehanna River, VA
        [-74.94442, 40.34478, np.nan, np.nan],  # Delaware River, NJ
        [-78.425288, 34.508177, np.nan, np.nan],  # Cape Fear River, NC
        [-91.72306, 31.04462, np.nan, np.nan],  # Red River (upstream of Atchafalaya River), LA
        [-80.10808, 33.50005, np.nan, np.nan],  # Santee River, SC
        [-79.81703, 33.59694, np.nan, np.nan],  # Black River, SC
        [-79.57210, 33.71223, np.nan, np.nan],  # Black Mingo Creek, SC
        [-79.49997, 33.84686, np.nan, np.nan],  # Lynches River, SC
        [-79.48467, 33.93939, np.nan, np.nan],  # Pee Dee River, SC
        [-79.33247, 33.98196, np.nan, np.nan],  # Little Pee Dee River, SC
        [-77.917829, 34.749979, np.nan, np.nan],  # Northeast Cape Fear River, NC
        [-87.9523, 30.8472, np.nan, np.nan],  # Mobile River, AL
        [-96.974415, 28.673913, -97.002529, 28.750011],  # Blue Bayou, TX
        [-96.695401, 28.968284, -96.69652166667, 28.990345],  # Lavaca River, TX
        [-96.548436, 28.999706, -96.554498, 29.024612666667],  # Lake Texana, TX
        [-93.83342666667, 30.355123333333, -93.83342666667, 30.355123333333],  # Cypress Creek, TX
        [-89.764476, 30.551926, -89.76781133333, 30.538070666667],  # Lotts Creek, LA
        [-87.219805, 30.567296, -87.24471466667, 30.601442333333],  # Escambia River, FL
        [-83.987035, 30.331327, np.nan, np.nan],  # Horsehead Creek and Little River, FL
        [-83.928038, 30.30404, np.nan, np.nan],  # Bailey Mill Creek, FL
        [-82.950913, 29.958097, -82.99605566667, 30.007415],  # Suwannee River, FL
        
    ]).reshape(-1, 4)

    relocate_map = None  # if not None, use the prepared relocate_map to relocate sources

    relocate_sources(
        old_ss_dir=old_ss_dir,
        feeder_info_file=feeder_info_file,
        hgrid_fname=hgrid_fname,
        outdir=outdir,
        max_search_radius=max_search_radius,
        mandatory_sources_coor=mandatory_sources_coor,
        relocate_map=relocate_map
    )