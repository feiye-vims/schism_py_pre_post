import numpy as np
from pylib import read
vg = read('/sciclone/home/feiye/vgrid.in')
vg2 = read('/sciclone/home/feiye/TEMP/vgrid.in')
assert np.equal(vg.kbp, vg2.kbp).all()
assert np.allclose(vg.sigma, vg2.sigma)
pass

import numpy as np

import netCDF4



def cal_disturbance(
    depth: np.ndarray, elev: np.ndarray,
    city_node_idx_file: str = None,
    fillvalue: float = -99999.0
):
    """
    Calculate the maximum disturbance

    - inputs:
        depth: np.ndarray of shape (npoints,)
        elev: np.ndarray of shape (ntimes, npoints)
    - outputs:
        maxdist: np.ndarray of shape (npoints,)
    """
    if fillvalue > -1000:
        raise ValueError("fillvalue should be a large negative number")

    elev[np.isnan(elev)] = fillvalue  # deprecated, adcirc format

    disturb = elev.copy()  # same as elev in the ocean, so initialize with elev

    # read handle city node indices
    if city_node_idx_file is not None:
        city_node_idx = np.loadtxt(city_node_idx_file, encoding='utf-8').astype(bool)
    else:
        city_node_idx = np.zeros_like(depth, dtype=bool)
    # define land nodes, including city nodes
    land_node_idx = (depth < 0) | city_node_idx

    # Define disturbance:
    # On land, disturbance is the sum of elevation and depth, i.e., the water depth.
    # Also, disturbance is zero if the water depth is negative.
    disturb[:, land_node_idx] = np.maximum(0, elev[:, land_node_idx] + depth[land_node_idx])

    max_disturb = np.max(disturb, axis=0)
    time_idx_max_disturb = np.argmax(disturb, axis=0)

    # mask small max disturbance (< 0.3 m) on land (including cities)
    small_dist_on_land = (max_disturb < 0.3) * land_node_idx  # True if both conditions are met
    max_disturb[small_dist_on_land] = fillvalue

    return max_disturb, time_idx_max_disturb, disturb


import os
import xarray as xr

ncfiles = ['./out2d_20240911.nc', './out2d_20240912.nc', './out2d_20240913.nc']
for ncfile in ncfiles:
    print(f'Processing {ncfile}')
    os.system(f'cp ./Backup/{ncfile} {ncfile}')
    ds = xr.open_dataset(ncfile, mode='r+')
    elev = np.array(ds.variables['elevation'])
    depth = np.array(ds.variables['depth'])

    # _, _, disturbance = cal_disturbance(
    #     depth=depth, elev=elev,
    #     city_node_idx_file='/sciclone/schism10/feiye/STOFS3D-v7/Shared_with_NOAA/v7/Shared_for_CERA/inputs/city_poly.node_id.shadow.txt'
    # )

    elev[elev + depth <= 0.3] = -99999.0
    ds['elevation'].loc[...] = elev
    ds.to_netcdf(ncfile, mode='a')

    ds.close()

file = '/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v23.3/stofs_3d_atl.t12z.vsource.th_VIMS7.1'
data = np.loadtxt(file)


import numpy as np

# Create sample data
working_dir = './'  # Change to your preferred directory
times = np.arange(3600, 3600 * 6, 3600)  # Example times: 3600, 7200, 10800, ...
sources = np.random.random((len(times), 1)) * 10  # Random values for sources

# Save using fmt='%10.4f' (with leading spaces)
np.savetxt(f'{working_dir}/vsource_with_leading_spaces.th',
           np.c_[times, sources[:, 0]], fmt='%10.4f', delimiter=' ')

# Save using fmt='%.4f' (without leading spaces)
np.savetxt(f'{working_dir}/vsource_without_leading_spaces.th',
           np.c_[times, sources[:, 0]], fmt='%.4f', delimiter=' ')

# Save using fmt='%d %.4f' for mixed formatting
np.savetxt(f'{working_dir}/vsource_mixed_formatting.th',
           np.c_[times, sources[:, 0]], fmt='%d %.4f', delimiter=' ')

print("Files saved:")
print("- vsource_with_leading_spaces.th")
print("- vsource_without_leading_spaces.th")
print("- vsource_mixed_formatting.th")


pass
