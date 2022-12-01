import numpy as np
from schism_py_pre_post.Grid.SourceSinkIn import source_sink
from glob import glob 
from scipy import spatial
import xarray as xr
from pylib import schism_grid
import pandas as pd


def nearest_neighbour(points_a, points_b):
    tree = spatial.cKDTree(points_b)
    return tree.query(points_a)[1]

def Sflux_var_at_points(sflux_dir, sflux_set='air_1', points=None, var_name=None, nt=None):
    '''
    extract variable values at specified points from SCHISM-format sflux files
    '''

    sflux_files = sorted(glob(f'{sflux_dir}/sflux_{sflux_set}.*.nc'))
    # read the first sflux file
    sflux = xr.open_dataset(sflux_files[0])
    if nt is None:
        nt = len(sflux['time'].data)
    sflux.close()

    idx = nearest_neighbour(np.c_[points[:, 0], points[:, 1]], np.c_[sflux.lon.data.ravel(), sflux.lat.data.ravel()])
    sflux_files = sorted(glob(f'{sflux_dir}/sflux_{sflux_set}.*.nc'))
    nt_all = nt * len(sflux_files)
    npts = len(idx)

    var_at_pts = np.empty((nt_all, npts), dtype=float) * np.nan
    time = np.empty((nt_all), dtype='datetime64[ns]')

    for i, sflux_file in enumerate(sflux_files):
        sflux = xr.open_dataset(sflux_file)
        var = sflux[var_name].data
        var = var.reshape(*var.shape[:-2], -1)
        var_at_pts[i*nt:(i+1)*nt, :] = var[:nt, idx]
        time[i*nt:(i+1)*nt] = sflux.time.data[:nt]
        print(f'processing {sflux_file}')
        sflux.close()

    return [time, var_at_pts]


if __name__ == "__main__":
    rundir = '/sciclone/scr10/feiye/RUN02/'
    model_start_time_str = "2018-08-01 00:00:00"

    # read hgrid and get element info
    gd = schism_grid(f'{rundir}/hgrid.ll')
    gd.compute_ctr()
    lon_ctr, lat_ctr = gd.xctr, gd.yctr

    gd.x, gd.y = gd.proj(prj0='epsg:4326', prj1='epsg:26918')
    gd.compute_area()

    # read the original source/sink files
    orignial_source_sink = source_sink(source_dir=rundir, start_time_str=model_start_time_str)

    # make another set of source/sink files based on sflux
    # read variable values from sflux
    var_at_ele_center = Sflux_var_at_points(f'{rundir}/sflux/', sflux_set='prc_1', points=np.c_[lon_ctr, lat_ctr], var_name='prate')
    # read time info from sflux
    times = pd.to_datetime(var_at_ele_center[0].astype(str))
    start_time_str = times[0].strftime("%Y-%m-%d %H:%M:%S")
    timedeltas = (times-times[0]).total_seconds().astype(int)
    nt = times.shape[0]
    # calculate new sources for all elements: prate x area
    sources_sflux = np.zeros((nt, gd.ne), dtype=float)
    for it in range(nt):  # for each time
        sources_sflux[it, :] = 1e-3 * var_at_ele_center[1][it] * gd.area  # prate is in kg/m^2/s, i.e., 1e-3 m3/m2/s
    # make new source sink 
    added_ss = source_sink(
        source_dir=None,
        source_eles=(np.array(range(gd.ne))+1).tolist(),
        sink_eles=[1],  # dummy
        start_time_str=start_time_str,
        timedeltas=timedeltas,
        vsource_data=sources_sflux
    )
    added_ss.update_vars()

    # add original source and sflux source
    total_ss = orignial_source_sink + added_ss
    if np.isnan(total_ss.vsource.data.astype(float)).any():
        raise Exception('nan found in sources')
    if total_ss.vsource.data.min() < 0.0:
        raise Exception('negative sources found in sources')

    if np.isnan(total_ss.vsink.data.astype(float)).any():
        raise Exception('nan found in sinks')
    if total_ss.vsink.data.max() > 0.0:
        raise Exception('positive sinks found in sources')
    
    # write
    total_ss.writer(f'{rundir}/Orignial+sflux_source_sink/')
    # total_ss.nc_writer(f'{rundir}/Orignial+sflux_source_sink/')  # needs to adapt to the latest pylib 

    # import pickle
    # with open('tmp.pkl', 'wb') as file:
    #     pickle.dump([orignial_source_sink, added_ss], file)
    # with open('tmp.pkl', 'rb') as file:
    #     [orignial_source_sink, added_ss] = pickle.load(file)

