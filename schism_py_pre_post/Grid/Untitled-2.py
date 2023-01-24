#%%
import xarray as xr



#%%
from pylib import schism_grid
from netCDF4 import Dataset
import matplotlib.pyplot as plt
# %%
gd = schism_grid('/sciclone/schism10/feiye/STOFS3D-v4/Inputs/v11.7/hgrid.ll.pkl')
# %%
schout = Dataset('/sciclone/home10/feiye/schout_adcirc_20220518.nc')
gd.dp = schout.variables['disturbance_max'][:]
pass
# %%
gd.plot(fmt=1, clim=[0.5, 2.8], cmap='jet')
plt.show()
# %%

