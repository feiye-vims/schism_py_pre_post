import xarray as xr
import numpy as np


ncfile = '/sciclone/schism10/feiye/STOFS3D-v7/Shared_with_NOAA/v7/Shared_for_CERA/extract/schout_adcirc_1.nc'
ncfile = '/sciclone/schism10/feiye/STOFS3D-v7/Shared_with_NOAA/v7/Shared_for_CERA/outputs/schout_adcirc_20240915.nc'
ncfile = '/sciclone/schism10/feiye/STOFS3D-v7/Shared_with_NOAA/v7/Shared_for_CERA/outputs/schout_adcirc_20240929.nc'

nc = xr.open_dataset(ncfile)
elev = np.array(nc['zeta'])
zeta_max = np.array(nc['zeta_max'])

id = 2498754
assert np.allclose(np.max(elev[:, id]), zeta_max[id]), 'zeta_max is not the max of zeta'

pass
