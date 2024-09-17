import matplotlib.pyplot as plt
from schism_py_pre_post.Timeseries.TimeHistory import TimeHistory
from pylib_experimental.schism_file import cread_schism_hgrid
from schism_py_pre_post.Grid.SourceSinkIn import SourceSinkIn, source_sink
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

    # vs = TimeHistory('/sciclone/schism10/feiye/STOFS3D-v4/Inputs/Iv4/20220502/vsource.th')
    # ss = SourceSinkIn('/sciclone/schism10/feiye/STOFS3D-v4/Inputs/Iv4/20220502/source_sink.in')
    # ss.toPropFile(ne=5654161)
    # vs1 = TimeHistory('/sciclone/schism10/feiye/STOFS3D-v4/RUN23h/vsource.th')
    # ss1 = SourceSinkIn('/sciclone/schism10/feiye/STOFS3D-v4/RUN23h/source_sink.in')
    # ms = TimeHistory('/sciclone/home10/feiye/ChesBay/RUN110y/msource.th')
    # ss1 = source_sink('/sciclone/schism10/feiye/Coastal_Act/Pumps/11LL_3days/')
    # vs = TimeHistory('/sciclone/schism10/feiye/STOFS3D-v4/fcst_run/vsource.th')

    WDIR = '/sciclone/schism10/feiye/STOFS3D-v8/I03f/Relocated_SS/'
    # ss = source_sink(WDIR)
    # ss.vsource.df['131606'].plot()

    [max_val, avg_source_sink, fig] = viz_source(WDIR, '2022-05-01 00:00:00', i_show_plot=-1, scale=1.e3)

    pass
