"""
Plot 2D variables (or a slab/layer of 3D variables) from SCHISM outputs,
use mpi to parallelize the plotting of multiple time stamps
"""
import os
import gc
import copy
from mpi4py import MPI
from datetime import datetime, timedelta

import numpy as np
import matplotlib
# from mpl_toolkits.basemap import Basemap  # mamba install basemap
import matplotlib.pyplot as plt
import netCDF4

import cartopy
import cartopy.geodesic as cgeo
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import io
from urllib.request import urlopen, Request
from PIL import Image

from pylib_experimental.schism_file import cread_schism_hgrid
from pylib import grd2sms, schism_grid


try:
    matplotlib.use('TkAgg')
except:
    matplotlib.use('Agg')

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
# print(f'process {rank} of {size}\n')


def my_mpi_idx(N, size, rank):
    '''
    Distribute N tasks to {size} ranks.
    The return value is a bool vector of the shape (N, ),
    with True indices indicating tasks for the current rank.
    '''

    # initialize a bool array of size N
    # for the current rank to decide which tasks to handle
    # (True for the tasks to handle, False for the tasks to skip)
    my_idx = np.zeros((N, ), dtype=bool)

    if N <= size:  # trivial case: more ranks than tasks
        if rank < N:
            my_idx[rank] = True

    else:  # cyclic distribution
        for i in range(N):
            if (i % size) == rank:
                my_idx[i] = True

    # my_idx = np.zeros((N, ), dtype=bool)
    # n_per_rank, _ = divmod(N, size)
    # n_per_rank = n_per_rank + 1
    # my_idx[rank*n_per_rank:min((rank+1)*n_per_rank, N)] = True

    return my_idx

def schism2sms_parallel(
    rundir='./',
    model_start_time = datetime.strptime('2014-12-01 00:00:00', "%Y-%m-%d %H:%M:%S"),
    var_str=None, var_proc=None, stacks=[1, 2, 3],
    output_dir=None, snapshots_times=(),
    iOverWrite=False
):
    if var_str in ['elevation']:
        fname_prefix = 'out2d'
    else:
        fname_prefix = var_str

    # devide tasks
    stacks_proc = stacks[my_mpi_idx(len(stacks), size, rank)]

    # stacks_proc = np.flip(stacks_proc)
    print(f'process {rank} handles {stacks_proc}\n')

    os.makedirs(output_dir, exist_ok=True)

    gd = cread_schism_hgrid(f'{rundir}/hgrid.gr3')

    for stack in stacks_proc:
        fname = f'{rundir}/outputs/{fname_prefix}_{stack}.nc'
        print(f'Core {rank} reading {fname}\n')
        my_nc = netCDF4.Dataset(fname)

        this_time_dts = my_nc.variables['time']

        for i_time, this_time_dt in enumerate(this_time_dts):
            this_time = model_start_time + timedelta(seconds=this_time_dt.data.item())
            # save a copy as *.2dm
            if this_time in snapshots_times:
                # make a copy of the grid because the grid will be modified
                gd_copy = copy.deepcopy(gd)
                if 'var' not in locals():
                    var = np.array(my_nc.variables[var_str])
                if var_proc is not None:
                    var, title_var_str = var_proc(gd_copy, var)
                else:
                    title_var_str = var_str

                title_str = f"{title_var_str} {this_time.strftime('%Y-%m-%d %H:%M:%S')}"
                savefilename = f'{output_dir}/{title_str.replace(" ", "_")}.2dm'
                if not iOverWrite and os.path.exists(savefilename):
                    print(f'Core {rank} skipping: {savefilename}\n')
                    continue  # skip existing figures
                if var.ndim == 3:
                    gd_copy.dp = var[i_time, :, -1]
                elif var.ndim == 2:
                    gd_copy.dp = var[i_time, :]
                    gd_copy.dp[np.isnan(gd_copy.dp)] = -9999
                grd2sms(gd_copy, savefilename)

        my_nc.close()
        if 'var' in locals():
            del var
        gc.collect()

    print(f'Core {rank} finishing ...\n\n\n')
    comm.Barrier()


def plot_schism2D_parallel(
    rundir='./',
    model_start_time=datetime.strptime('2014-12-01 00:00:00', "%Y-%m-%d %H:%M:%S"),
    var_str=None, var_proc=None, stacks=None, time_steps=None,
    plot_params: dict = None, output_dir=None, iOverWrite=False
):
    """
    function for plotting 2D variables of the schism outputs;
    mpi can be used when plotting many time stamps
    """
    if plot_params is None:
        plot_params = {
            'xlim': [-92.3, -88.5], 'ylim': [28.8, 31.2], 'clim': [-2, 10],
        }
    if stacks is None:
        stacks = [1, 2]
    if time_steps == []:
        time_steps = None
    if var_str in ['elevation']:
        fname_prefix = 'out2d'
    else:
        fname_prefix = var_str

    # devide tasks
    stacks_proc = stacks[my_mpi_idx(len(stacks), size, rank)]

    # stacks_proc = np.flip(stacks_proc)
    print(f'process {rank} handles {stacks_proc}\n')

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    gd = cread_schism_hgrid(f'{rundir}/hgrid.gr3')

    for stack in stacks_proc:
        fname = f'{rundir}/outputs/{fname_prefix}_{stack}.nc'
        print(f'Core {rank} reading {fname}\n')
        my_nc = netCDF4.Dataset(fname)
        var = np.array(my_nc.variables[var_str])

        gd_copy = copy.deepcopy(gd)  # get a copy of the original grid because the copy will be modified
        if var_proc is not None:
            var, var_title_str = var_proc(gd_copy, var)

        this_time_dts = my_nc.variables['time']
        if time_steps is None:
            time_steps = range(len(this_time_dts))  # all time stamps

        for j, this_time_dt in enumerate(this_time_dts):
            if j not in time_steps:  # only plot selected time stamps
                continue
            this_time = model_start_time + timedelta(seconds=this_time_dt.data.item())
            title_str = f"{var_title_str} {this_time.strftime('%Y-%m-%d %H:%M:%S')}"
            savefilename = f'{output_dir}/{title_str.replace(" ", "_")}.png'

            if not iOverWrite and os.path.exists(savefilename):
                print(f'Core {rank} skipping: {savefilename}\n')
                continue  # skip existing figures

            print(f'Core {rank} plotting time: {title_str} from {fname}\n\n')
            if var.ndim == 3:
                gd_copy.dp = var[j, :, -1]
            elif var.ndim == 2:
                gd_copy.dp = var[j, :]

            crs = ccrs.PlateCarree()
            _, ax = plt.subplots(figsize=(12, 6), subplot_kw={'projection': crs})

            extent = [plot_params['xlim'][0], plot_params['xlim'][1], plot_params['ylim'][0], plot_params['ylim'][1]]
            ax.set_extent(extent, crs=crs)

            ibasemap = False
            if ibasemap:
                cimgt.OSM.get_image = image_spoof # reformat web request for street map spoofing
                img = cimgt.OSM() # spoofed, downloaded street map

                cimgt.QuadtreeTiles.get_image = image_spoof
                img = cimgt.QuadtreeTiles()  # satelite image

                ax.add_image(img, 14)
                ax.coastlines()
                ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

                # m = Basemap(projection="merc", resolution="f")  # mamba install basemap-data-hires
                # m.shadedrelief()
                # m.drawcoastlines()

            gd_copy.plot_grid(ax=ax, fmt=1, levels=21, cmap='jet', **plot_params)
            if ibasemap:
                gl = ax.gridlines(draw_labels=True, crs=crs,
                        color='k',lw=0.5)
                gl.top_labels = False
                gl.right_labels = False
                gl.xformatter = cartopy.mpl.gridliner.LONGITUDE_FORMATTER
                gl.yformatter = cartopy.mpl.gridliner.LATITUDE_FORMATTER

            plt.title(title_str)
            plt.gca().set_aspect('equal', 'box')
            plt.savefig(savefilename, dpi=400)
            # plt.show()
            plt.clf()
            plt.close()

        my_nc.close()
        del var, my_nc
        gc.collect()

    print(f'Core {rank} finishing ...\n\n\n')
    comm.Barrier()

def image_spoof(self, tile):
    '''this function reformats web requests from OSM for cartopy
    Heavily based on code by Joshua Hrisko at:
        https://makersportal.com/blog/2020/4/24/geographic-visualizations-in-python-with-cartopy'''

    url = self._image_url(tile)                # get the url of the street map API
    req = Request(url)                         # start request
    req.add_header('User-agent','Anaconda 3')  # add user agent to request
    fh = urlopen(req)
    im_data = io.BytesIO(fh.read())            # get image
    fh.close()                                 # close url
    img = Image.open(im_data)                  # open image with PIL
    img = img.convert(self.desired_tile_form)  # set image format
    return img, self.tileextent(tile), 'lower' # reformat for cartopy


def mask_dry_nodes(hgrid_obj, var, var_str, elevation, min_depth=0.01):
    '''mask dry nodes in the variable
    min_depth: minimum depth to be considered as water
    '''
    var = np.array(var)
    elev = np.array(elevation)

    for i_time in range(var.shape[0]):
        dry_idx = hgrid_obj.dp + elev[i_time, :] < min_depth
        var[i_time, dry_idx] = np.nan

    return var, 'masked_' + var_str

def mask_dry_elevation(hgrid_obj, elevation):
    return mask_dry_nodes(hgrid_obj, elevation, 'elevation', elevation, min_depth=0)

def get_disturbance(hgrid_obj, elevation):
    '''elevation can have a time dimension,
    '''
    elev = np.array(elevation)
    disturbance = np.array(elevation)  # elevation for water nodes

    land = hgrid_obj.dp < 0.0
    for i_time in range(disturbance.shape[0]):
        disturbance[i_time, land] = elev[i_time, land] + hgrid_obj.dp[land]  # inundation depth for land nodes

    return disturbance, 'disturbance'

def get_positive_disturbance(hgrid_obj, elevation):
    disturbance, _ = get_disturbance(hgrid_obj, elevation)

    disturbance[disturbance<0] = np.nan

    return disturbance, 'positive_disturbance'

plot_param_dict = {
    'Pearl River': {
        'xlim': [-89.87, -89.55], 'ylim': [30.30, 30.55], 'clim': [-2, 10],
    },
    'LA': {
        'xlim': [-92.3, -88.5], 'ylim': [28.8, 31.2], 'clim': [-2, 10],
    },
    'New Orleans': {
        'xlim': [-90.19665, -89.96121], 'ylim': [29.89496, 30.05080], 'clim': [-2, 10],
    },
    'Outfall Canal': {
        'xlim': [-90.10806, -90.00778], 'ylim': [29.97116, 30.03879], 'clim': [-2, 10],
    }
}

if __name__ == "__main__":
    # sample inputs
    RUNDIR = '/sciclone/schism10/feiye/STOFS3D-v8/R10a/'
    output_dir = f'{RUNDIR}/outputs/'
    model_start_time = datetime.strptime('2024-03-05', "%Y-%m-%d")
    VAR_STR = 'elevation'
    var_proc = mask_dry_elevation
    stacks = np.arange(1, 36)  # must cover the plot time stamps
    time_steps = []  # None (all time steps), or a list of time steps to plot
    plot_params = plot_param_dict['Pearl River']

    # -------------------- sample outputing to *.2dm -------------------------
    snapshots_times = (
        datetime.strptime('2024-03-06 00:00:00', "%Y-%m-%d %H:%M:%S"),
        datetime.strptime('2024-03-13 00:00:00', "%Y-%m-%d %H:%M:%S"),
        datetime.strptime('2024-03-20 00:00:00', "%Y-%m-%d %H:%M:%S"),
        datetime.strptime('2024-03-27 00:00:00', "%Y-%m-%d %H:%M:%S"),
    )
    if snapshots_times:
        schism2sms_parallel(
            RUNDIR, model_start_time, VAR_STR, var_proc, stacks,
            output_dir, snapshots_times, iOverWrite=True)

    # ----------------------------- sample generating plot -------------------------
    plot_schism2D_parallel(
        RUNDIR, model_start_time, VAR_STR, var_proc, stacks,
        time_steps, plot_params, output_dir, iOverWrite=False)

    print('All done!')
