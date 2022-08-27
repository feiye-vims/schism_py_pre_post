# %%
from pylib import schism_grid, read_schism_vgrid
import os
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import pickle
from IPython.display import set_matplotlib_formats
set_matplotlib_formats('svg')


# %%
gd_fname = '/sciclone/schism10/feiye/STOFS3D-v4/Inputs/I23o/hgrid.gr3'
gd_cache_fname = os.path.splitext(gd_fname)[0] + '.pkl'
if os.path.exists(gd_cache_fname):
    gd = schism_grid(gd_cache_fname)
else:
    gd = schism_grid(gd_fname)
    gd.save(gd_cache_fname)

# %%
vg_fname = '/sciclone/schism10/feiye/STOFS3D-v4/RUN23o/vgrid.in'
vg_cache_fname = os.path.splitext(vg_fname)[0] + '.pkl'
if os.path.exists(vg_cache_fname):
    with open(vg_cache_fname, 'rb') as handle:
        vg = pickle.load(handle)
else:
    vg = read_schism_vgrid(vg_fname)
    with open(vg_cache_fname, 'wb') as handle:
        pickle.dump(vg, handle, protocol=pickle.HIGHEST_PROTOCOL)


# %%
var_name = "salinity"
it = -1
isurf = False
caxis = [0, 35]; xlim = [-77, -75]; ylim = [37, 39]
plt.figure(figsize=(7, 7))

fname = f'/sciclone/schism10/feiye/STOFS3D-v4/RUN23p9/outputs/{var_name}_78.nc'
fname = f'/sciclone/scr10/lcui01/ICOGS3D/outputs_RUN23p1/{var_name}_21.nc'

my_nc = netCDF4.Dataset(fname)
var = my_nc.variables[var_name]

if isurf:
    value = var[it, :, -1]
else:
    value = var[it, :, :]
    value = value[np.arange(gd.np), vg.kbp]

gd.plot_grid(fmt=1, value=value, clim=caxis, levels=31, ticks=np.arange(caxis[0], caxis[1], 2), cmap='jet', xlim=xlim, ylim=ylim)
plt.show()
my_nc.close()

# %%
pass
