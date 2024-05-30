#%%
import os
from glob import glob
from pathlib import Path

from matplotlib import pyplot as plt
import netCDF4
import xarray as xr
import numpy as np

from pylib_experimental.schism_file import cread_schism_hgrid
from tqdm import tqdm


#%%

# compare two nc files
nc1 = xr.open_dataset('/sciclone/schism10/feiye/STOFS3D-v7/Shared_with_NOAA/v7/Shared_for_CERA/extract/1.nc')
nc2 = xr.open_dataset('/sciclone/schism10/feiye/STOFS3D-v7/Shared_with_NOAA/v7/Shared_for_CERA/extract/schout_UV4.5m_1.nc')

if sorted(list(nc1.variables.keys())) != sorted(list(nc2.variables.keys())):
    raise ValueError('Variables are not the same')
var_list = list(nc1.variables.keys())

for key in var_list:
    dtype = nc1[key].dtype
    print(f'Comparing {key}, dtype: {dtype}')
    if np.issubdtype(dtype, np.floating):
        is_equal = np.allclose(np.array(nc1[key]), np.array(nc2[key]), equal_nan=True, atol=1e-7, rtol=1e-7)
    else:
        is_equal = nc1[key].equals(nc2[key])
    if not is_equal:
        raise ValueError(f'{key} is not equal')

    # close_idx = np.isclose(np.array(nc1[key]), np.array(nc2[key]), equal_nan=True, atol=1e-6, rtol=1e-6)

#%%
# read from local
hgrid_obj = cread_schism_hgrid('/sciclone/schism10/feiye/STOFS3D-v7/Inputs/I14b/hgrid.gr3')
with xr.open_dataset('/sciclone/home/feiye/s1/STOFS3D-v7/Shared_with_NOAA/v7/Shared_for_CERA/extract/schout_UV4.5m_1.nc') as ds:
    u = np.array(ds['uvel4.5'])
    v = np.array(ds['vvel4.5'])
    var = np.sqrt(u**2 + v**2)
    # var = np.array(ds['disturbance_max'])

hgrid_obj.plot(value=var[12, :], fmt=1, clim=[0, 2], cmap='jet')
plt.title('Velocity magnitude near bottom')
plt.show()
pass

#%%
# read from aws
import s3fs
s3_path = 's3://noaa-nos-stofs3d-pds/STOFS-3D-Atl/stofs_3d_atl.20240325/stofs_3d_atl.t12z.f001_024.field2d.nc'
fs = s3fs.S3FileSystem(anon=True)  # Set anon to True if the bucket is public
with fs.open(s3_path, mode='rb') as f:
    with netCDF4.Dataset('dummy', mode='r', memory=f.read()) as ds:
        print(ds.variables)  # Print list of variables
        # Access data from a variable, for example:
        variable_data = ds.variables['salt_surface'][:]
        print(variable_data)

pass

#%%
# change the precipitation in all nc files
files = sorted(glob('/sciclone/schism10/feiye/STOFS3D-v7/Inputs/I17f/Sflux/I17a_sflux/sflux_prc_2.*nc'))
outdir = '/sciclone/schism10/feiye/STOFS3D-v7/Inputs/I17f/Sflux/sflux/'
for file in tqdm(files):
    fname = Path(file).name
    os.system(f'cp {file} {outdir}/{fname}')
    with netCDF4.Dataset(f'{outdir}/{fname}', 'r+') as ds:
        ds['prate'][:] = ds['prate'][:] * 0.516

pass

#%%
# change the base year in all nc files
original_year = '2005'
new_year = '2021'
outdir = '/sciclone/schism10/feiye/STOFS3D-v7/Inputs/I17a/Sflux/sflux/'
files = glob('/sciclone/schism10/feiye/STOFS3D-v7/Inputs/I17a/Sflux/RUN20_sflux/sflux_???_2.*nc')
for file in tqdm(files):
    fname = Path(file).name
    os.system(f'cp {file} {outdir}/{fname}')
    with netCDF4.Dataset(f'{outdir}/{fname}', 'r+') as ds:
        based_date = ds['time'].base_date
        based_date[0] = int(new_year)
        ds['time'].base_date = based_date
        ds['time'].units = ds['time'].units.replace(original_year, new_year)
        
# plot velocity vector at a point
u_files = sorted(glob('/sciclone/schism10/feiye/STOFS3D-v7/Runs/R15/outputs/horizontalVelX_*.nc'))
v_files = sorted(glob('/sciclone/schism10/feiye/STOFS3D-v7/Runs/R15/outputs/horizontalVelY_*.nc'))

poi = [-76.52, 39.22]
hgrid_obj = cread_schism_hgrid('/sciclone/schism10/feiye/STOFS3D-v7/Runs/R15/hgrid.gr3')
idx = np.argmin((hgrid_obj.x - poi[0])**2 + (hgrid_obj.y - poi[1])**2)

for i, [u_file, v_file] in tqdm(enumerate(zip(u_files, v_files)), total=len(u_files)):
    with xr.open_dataset(u_file) as u_data:
        u_poi = np.array(u_data['horizontalVelX'])[:, idx, -1]
    with xr.open_dataset(v_file) as v_data:
        v_poi = np.array(v_data['horizontalVelY'])[:, idx, -1]

    if i == 0:
        time = np.array(u_data['time'])
        # initialize arrays
        u = np.zeros(len(time)*len(u_files), )
        v = np.zeros(len(time)*len(u_files), )
        t = np.zeros((len(time)*len(u_files), ), dtype='datetime64[s]')
        
    u[i*len(time):(i+1)*len(time)] = u_poi
    v[i*len(time):(i+1)*len(time)] = v_poi
    t[i*len(time):(i+1)*len(time)] = time

plt.quiver(t, np.zeros_like(t), u, v, width = 0.003, scale=0.001, angles='xy', scale_units='xy')
plt.show()