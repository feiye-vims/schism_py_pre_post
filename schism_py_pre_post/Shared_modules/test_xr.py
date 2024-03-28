import xarray as xr
import numpy as np

from glob import glob

files = glob('/sciclone/schism10/feiye/STOFS3D-v7/Runs/R14c/sflux/sflux_air_1.*.nc')
nc = xr.open_mfdataset(files, combine='nested', concat_dim='time')
print(np.max(abs(nc['uwind'].values + 1j* nc['vwind'].values)))


nc1 = xr.open_dataset('/sciclone/schism10/feiye/STOFS3D-v7/Inputs/I14c/Bnd/elev2D.th.nc')
nc2 = xr.open_dataset('/sciclone/schism10/lcui01/schism20/ICOGS/ICOGS3D/Scripts/RUN13u/bnd/elev2D.th.nc.aviso')

var1 = nc1['time_series'].values
var2 = nc2['time_series'].values

var3 = var2.copy()
var3[100,100,0,0] = 0
np.isclose(var1, var3, atol=1e-8).all()