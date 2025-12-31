import matplotlib.pyplot as plt
from schism_py_pre_post.Timeseries.TimeHistory import TimeHistory
from pylib import read
from pylib_experimental.schism_file import cread_schism_hgrid, source_sink, SourceSinkIn
import numpy as np
import os


def find_duplicates(arr):
    """
    Find duplicate values and their indices in a 1D numpy array.
    
    Parameters:
    arr (numpy array): Input 1D numpy array.
    
    Returns:
    dict: A dictionary where keys are duplicate values and values are lists of their indices.
    """
    # Dictionary to keep track of indices for each value
    index_dict = {}
    
    # Iterate over the array and store indices
    for idx, value in enumerate(arr):
        if value in index_dict:
            index_dict[value].append(idx)
        else:
            index_dict[value] = [idx]
    
    # Filter out entries with only one occurrence
    duplicates = {key: indices for key, indices in index_dict.items() if len(indices) > 1}
    
    return duplicates


def find_source_ele_by_coords(hgrid, source_sink_in, x, y):
    """
    Find nearest element index in the hgrid given x, y coordinates.
    x, y: float or array-like
    """
    from scipy.spatial import cKDTree
    hgrid.compute_ctr()
    source_x = hgrid.xctr[source_sink_in.ele_groups[0] - 1]
    source_y = hgrid.yctr[source_sink_in.ele_groups[0] - 1]
    tree = cKDTree(np.c_[source_x, source_y])
    dist, idx = tree.query(np.c_[x, y], k=1)
    return dist, idx


def plot_source_at_coords(hgrid, source_sink_in, vsource_th, x_coords, y_coords, tol=5e-2):
    '''
    Plot source/sink time series at given coordinates
    x_coords, y_coords: list of floats
    tol: tolerance to find nearest element; default is 5e-2 (about 5 km)
    '''
    dist, ele_idx = find_source_ele_by_coords(hgrid, source_sink_in, x_coords, y_coords)
    for i, idx in enumerate(ele_idx):
        if dist[i] > tol:
            print(f"Warning: No source/sink found within tolerance {tol} for point ({x_coords[i]}, {y_coords[i]})")
        source_ele_id = str(source_sink_in.ele_groups[0][idx])
        plt.plot(vsource_th.df.index, vsource_th.df[str(source_ele_id)], label=f'Point {i} at ({x_coords[i]}, {y_coords[i]})')

    plt.xlabel('Time')
    plt.ylabel('Source/Sink Value')
    plt.legend()
    plt.show()
    

def viz_source(w_dir, start_time_str, i_show_plot=0, scale=1e3, i_nc=False):
    '''
    Visualize sources/sinks on hgrid
        i_show_plot = 0: save figure but don't show; 1: show figure; -1: save *.xyz to be visualized with the preset symbols "sink_symbols.qgis" in Qgis (faster with large meshes)
    '''
    # initialize hgrid
    my_hgrid = cread_schism_hgrid(w_dir + '/hgrid.gr3')

    # draw hgrid bnd
    my_hgrid.plot()

    # read source/sink
    if i_nc:
        source_sink_fname = ['source.nc']
    else:
        source_sink_fname = ['vsource.th', 'vsink.th']
        source_sink_color = ['r', 'b']
    max_val = 0.0
    avg_source_sink = []

    for i in range(0, 2):  # 1st: source; 2nd: sink
        if os.path.exists(f'{w_dir}/{source_sink_fname[i]}'):
            my_ele_idx = SourceSinkIn(w_dir + '/source_sink.in', 2).ip_group[i] - 1
            my_hgrid.compute_ctr()
            my_ele_xyz = np.c_[my_hgrid.xctr[my_ele_idx], my_hgrid.yctr[my_ele_idx], my_hgrid.dpe[my_ele_idx]]

            # my_ele_th = TimeHistory(
            #     w_dir + source_sink_fname[i], start_time_str, -9999)
            # my_ele_xyz[:, 2] = np.mean(my_ele_th.data, axis=0)

            my_ele_th_data = np.loadtxt(f"{w_dir}/{source_sink_fname[i]}")
            my_ele_xyz[:, 2] = np.mean(my_ele_th_data[:, 1:], axis=0)

            print(f'max_val = {np.max(abs(my_ele_xyz[:, 2]))}')
            max_val = max(max_val, np.max(abs(my_ele_xyz[:, 2])))

            print(len(my_ele_idx))
            print(len(np.unique(my_ele_idx)))
            large = np.where(abs(my_ele_xyz[:, 2]) > 10)[0]
            print(len(large))
            print(len(np.unique(np.floor(my_ele_xyz[large, 2] * 1000))))

            fig = None
            if i_show_plot == -1:
                # put large values at the end to be visualized on top
                sorted_indices = np.argsort(my_ele_xyz[:, -1])
                my_ele_xyz = my_ele_xyz[sorted_indices]
                np.savetxt(f'{w_dir}/{source_sink_fname[i]}.xyz', my_ele_xyz)
            else:
                if scale > 0:
                    fig = plt.scatter(my_ele_xyz[:, 0], my_ele_xyz[:, 1],
                                      c=source_sink_color[i],
                                      s=abs(np.maximum(0.01, my_ele_xyz[:, 2]/max_val)*scale),
                                      alpha=0.7)
                else:
                    fig = plt.scatter(my_ele_xyz[:, 0], my_ele_xyz[:, 1],
                                      c=source_sink_color[i],
                                      s=abs(my_ele_xyz[:, 2]*0+10),  # neglecting volumes, only plotting locations, '10' is symbol size
                                      alpha=0.7)
                avg_source_sink.append(my_ele_xyz)
                my_ele_xyz = None

    if i_show_plot >= 0:
        plt.axis('equal')
        if i_show_plot == 1:
            plt.show()
        plt.savefig('1.png', dpi=600)
        if i_show_plot == 1:
            return [max_val, avg_source_sink, None]
        elif i_show_plot == 0:
            return [max_val, avg_source_sink, fig]
    else:
        return [max_val, avg_source_sink, None]


if __name__ == "__main__":
    '''Sample usage'''
    hgrid = read('/sciclone/schism10/feiye/STOFS3D-v7.3/R21l/hgrid.gr3')
    ss = source_sink.from_files('/sciclone/schism10/feiye/STOFS3D-v7.3/R21j/')
    xy = np.array([
        [-91.7223, 31.0457],  # Atchafalaya River
        [-91.56, 31.05],  # Mississippi River
        [-76.1096, 39.5874],  # Susquehanna River
        [-69.77256, 44.31494],  # Kennebec River, ME
        [-73.90869, 42.13509],  # Hudson River, NY
        [-74.94442, 40.34478],  # Delaware River, NJ
        [-78.425288, 34.508177],  # Cape Fear River, NC
        [-80.10808, 33.50005],  # Santee River, SC
        [-79.81703, 33.59694],  # Black River, SC
        [-79.57210, 33.71223],  # Black Mingo Creek, SC
        [-79.49997, 33.84686],  # Lynches River, SC
        [-79.48467, 33.93939],  # Pee Dee River, SC
        [-79.33247, 33.98196],  # Little Pee Dee River, SC
        [-77.917829, 34.749979],  # Northeast Cape Fear River, NC
        [-87.9523, 30.8472],  # Mobile River, AL
        [-96.695401, 28.968284],  # Lavaca River, TX
        [-96.548436, 28.999706],  # Lake Texana, TX
        [-93.83342666667, 30.355123333333],  # Cypress Creek, TX
        [-89.764476, 30.551926],  # Lotts Creek, LA
        [-87.219805, 30.567296],  # Escambia River, FL
        [-83.987035, 30.331327],  # Horsehead Creek and Little River, FL
        [-83.928038, 30.30404],  # Bailey Mill Creek, FL
        [-82.950913, 29.958097],  # Suwannee River, FL
        [-81.02370433333333, 27.315079666666666],  # Kissimmee River, FL
        [-81.997572, 30.786870],  # St Marys River, FL
        [-79.43425, 33.84487],  # Lyches River, SC
        [-74.74868, 39.47915],  # Great Egg Harbor River, NJ
        [-73.94009733333333, 42.06972966666667],  # Saugeties Creek, NY
        [-73.971293, 41.920595999999996],  # Hudson River branch, NY
        [-73.92918633333333, 41.592421333333334],  # Hudson River branch, NY
        [-73.07229533333333, 41.303546000000004],  # Housatonic River, CT
        [-72.625735, 41.656137666666666],  # Connecticut River, CT
        [-72.64970633333333, 41.572111666666665],  # Mattabesset River, CT
        [-72.470818, 41.47020933333334],  # Salmon River, CT
        [-72.11158266666666, 41.455657333333335],  # Stony Brook, CT
        [-72.090553, 41.535118000000004],  # Yantic River, CT
        [-72.06195833333334, 41.525600000000004],  # Quinebaug River, CT
    ])
    plot_source_at_coords(hgrid, ss.source_sink_in, ss.vsource, xy[:, 0], xy[:, 1], tol=5e-2)

    # vs = TimeHistory('/sciclone/schism10/feiye/STOFS3D-v4/Inputs/Iv4/20220502/vsource.th')
    # ss = SourceSinkIn('/sciclone/schism10/feiye/STOFS3D-v4/Inputs/Iv4/20220502/source_sink.in')
    # ss.toPropFile(ne=5654161)
    # vs1 = TimeHistory('/sciclone/schism10/feiye/STOFS3D-v4/RUN23h/vsource.th')
    # ss1 = SourceSinkIn('/sciclone/schism10/feiye/STOFS3D-v4/RUN23h/source_sink.in')
    # ms = TimeHistory('/sciclone/home10/feiye/ChesBay/RUN110y/msource.th')
    # ss1 = source_sink('/sciclone/schism10/feiye/Coastal_Act/Pumps/11LL_3days/')
    # vs = TimeHistory('/sciclone/schism10/feiye/STOFS3D-v4/fcst_run/vsource.th')

    WDIR = '/sciclone/schism10/feiye/STOFS3D-v8/I13x_v7/Source_sink/'
    ss = source_sink(WDIR)
    ss.vsource.df['5773144'].plot()
    plt.show()

    [max_val, avg_source_sink, fig] = viz_source(WDIR, '2024-03-05 00:00:00', i_show_plot=-1, scale=1.e3)

    pass
