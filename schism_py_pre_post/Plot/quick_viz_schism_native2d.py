from pylib_experimental.schism_file import cread_schism_hgrid
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from IPython.display import set_matplotlib_formats
from pathlib import Path
import xarray as xr
import copy
set_matplotlib_formats('svg')


gd_fname = '/sciclone/schism10/feiye/STOFS3D-v7/Shared_with_NOAA/v7/Shared_for_CERA/From_Zizang/hgrid.gr3'
gd = cread_schism_hgrid(gd_fname)

var_name = "zeta_max"  # "salt_surface"
fnames = ['/sciclone/schism10/feiye/STOFS3D-v7/Shared_with_NOAA/v7/Shared_for_CERA/From_Zizang/schout_adcirc_20250203.nc']
my_nc = xr.open_mfdataset(fnames)
value = np.array(my_nc[var_name])
caxis = [-2, 2]
plt.figure(figsize=(7, 7))
xlim = None
ylim=None

gd.plot(fmt=1, value=value, clim=caxis, levels=31, ticks=np.arange(caxis[0], caxis[1], 2), cmap='jet', xlim=xlim, ylim=ylim)
plt.gca().set_aspect('equal', 'box')
plt.show()
my_nc.close()

pass
