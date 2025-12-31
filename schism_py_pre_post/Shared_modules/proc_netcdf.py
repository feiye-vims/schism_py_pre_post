import os
import xarray as xr
import numpy as np
from copy import deepcopy
from pathlib import Path
from glob import glob

wdir = '/sciclone/schism10/feiye/STOFS3D-v8/I29k/InterpElev2D/'
for ncfile in ['uv3D.th.nc', 'SAL_3D.th.nc', 'TEM_3D.th.nc', 'elev2D.th.nc']:
    ds = xr.open_dataset(f'{wdir}/{ncfile}')

    time = ds['time'].values
    time_step = ds['time_step'].values[0]

    # realign time to start from zero and increase by time_step
    time = np.arange(len(time)) * time_step
    time0 = deepcopy(time)
    time0[1:] = time0[1:] - time[1] + time_step  # time[0] is always zero, so align to time[1], then add a time step
    assert np.allclose(time, time0), "Time realignment error"

    # update time in the original dataset
    ds = ds.assign_coords(time=('time', time))
    ds.to_netcdf(f'{wdir}/{Path(ncfile).stem}_realigned.nc')

    pass


output_dir = '/sciclone/home/feiye/s1/STOFS3D-v7.3/R21g/outputs'
ncfile_types = ['out2d', 'salinity', 'temperature', 'horizontalVelX', 'horizontalVelY', 'zCoordinates']
for ncfile_type in ncfile_types:
    ncfiles = sorted(glob(f'{output_dir}/{ncfile_type}_*.nc'))
    for i, ncfile in enumerate(ncfiles):
        os.symlink(ncfile, f'{wdir}/{ncfile_type}_{i+1}.nc')