import matplotlib.pyplot as plt
from schism_py_pre_post.Timeseries.TimeHistory import TimeHistory
from schism_py_pre_post.Grid.Hgrid_extended import read_schism_hgrid_cached
from schism_py_pre_post.Grid.SourceSinkIn import SourceSinkIn, source_sink
import numpy as np
import os


def viz_source(w_dir, start_time_str, i_show_plot=0, scale=1e3, i_nc=False):
    '''
    Visualize sources/sinks on hgrid
        i_show_plot = 0: save figure but don't show; 1: show figure; -1: save *.xyz to be visualized with the preset symbols "sink_symbols.qgis" in Qgis (faster with large meshes)
    '''
    # initialize hgrid
    my_hgrid = read_schism_hgrid_cached(w_dir + '/hgrid.gr3')

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
        if os.path.exists(w_dir + source_sink_fname[i]):
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

            fig = None
            if i_show_plot == -1:
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
        return [max_val, avg_source_sink, fig]
    else:
        return [max_val, avg_source_sink, None]


if __name__ == "__main__":
    # vs = TimeHistory('/sciclone/schism10/feiye/STOFS3D-v4/Inputs/Iv4/20220502/vsource.th')
    # ss = SourceSinkIn('/sciclone/schism10/feiye/STOFS3D-v4/Inputs/Iv4/20220502/source_sink.in')
    # vs1 = TimeHistory('/sciclone/schism10/feiye/STOFS3D-v4/RUN23h/vsource.th')
    # ss1 = SourceSinkIn('/sciclone/schism10/feiye/STOFS3D-v4/RUN23h/source_sink.in')
    # ms = TimeHistory('/sciclone/home10/feiye/ChesBay/RUN110y/msource.th')
    # ss1 = source_sink('/sciclone/schism10/feiye/Coastal_Act/Pumps/11LL_3days/')
    # vs = TimeHistory('/sciclone/schism10/feiye/STOFS3D-v4/fcst_run/vsource.th')
    '''
    ------viz source------
    '''
    w_dir = '/sciclone/schism10/feiye/STOFS3D-v6/Inputs/I15/Source_sink/Total_ss/'

    # ss = source_sink(run_dir)
    # ss.toPropFile(ne=5654161)

    [max_val, avg_source_sink, fig] = viz_source(w_dir, '2022-05-01 00:00:00', i_show_plot=-1, scale=1.e3)

    pass
