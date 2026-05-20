import xarray as xr
import numpy as np
from copy import deepcopy


wdir = '/sciclone/schism10/feiye/STOFS3D-v8/I35/InterpElev2D/'
for ncfile in ['SAL_3D.th.nc', 'TEM_3D.th.nc', 'uv3D.th.nc', 'elev2D.th.nc']:
    ds = xr.open_dataset(f'{wdir}/{ncfile}')

    time = ds['time'].values
    time_step = ds['time_step'].values[0]

    # realign time to start from zero and increase by time_step
    time = np.arange(len(time)) * time_step
    time0 = deepcopy(time)

    # interpolate_variables8.f90 duplicates the first time step and put it at time[0],
    # if we read from stack N, the first few time steps are 0, T, T+dt, T+2*dt, ...,
    # where T is the time of the first time step in the original dataset, and dt is the time step.
    # For example: 0, 11494800, 11498400, 11502000, 11505600, 11509200, 11512800, ...,
    # so we need to realign time to start from zero and increase by time_step as:
    # 0, 3600, 7200, 10800, 14400, 18000, 21600, ...
    time0[1:] = time0[1:] - time[1] + time_step  # time[0] is always zero, so align time[1:], then add a time step

    assert np.allclose(time, time0), "Time realignment error"

    # update time in the dataset and save to a new netcdf file
    ds = ds.assign_coords(time=('time', time))

    ds.to_netcdf(f'{wdir}/{ncfile.replace(".th.nc", "_realign.nc")}')
