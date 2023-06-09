# %%
from schism_py_pre_post.Grid.Hgrid_extended import read_schism_hgrid_cached
from pylib import read_schism_vgrid
import os
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import pickle
from IPython.display import set_matplotlib_formats
set_matplotlib_formats('svg')


# %%
gd_fname = '/sciclone/schism10/feiye/STOFS3D-v6/Inputs/I13u6/hgrid.gr3'
gd = read_schism_hgrid_cached(gd_fname)

# %%
vg_fname = '/sciclone/schism10/feiye/STOFS3D-v6/Inputs/I13u6/vgrid.in'
vg_cache_fname = os.path.splitext(vg_fname)[0] + '.pkl'
if os.path.exists(vg_cache_fname):
    with open(vg_cache_fname, 'rb') as handle:
        vg = pickle.load(handle)
else:
    vg = read_schism_vgrid(vg_fname)
    with open(vg_cache_fname, 'wb') as handle:
        pickle.dump(vg, handle, protocol=pickle.HIGHEST_PROTOCOL)


# %%
var_name = "temperature"
fname = f'/sciclone/home/feiye/Sync/{var_name}_275.nc'

it = -1
isurf = True
caxis = [0, 35]
xlim = None  # [-77, -75]
ylim = None  # [37, 40]
plt.figure(figsize=(7, 7))

my_nc = netCDF4.Dataset(fname)
var = my_nc.variables[var_name]

if isurf:
    value = var[it, :, -1]
else:
    value = var[it, :, :]
    value = value[np.arange(gd.np), vg.kbp]

gd.plot_grid(fmt=1, value=value, clim=caxis, levels=31, ticks=np.arange(caxis[0], caxis[1], 2), cmap='jet', xlim=xlim, ylim=ylim)
plt.gca().set_aspect('equal', 'box')
plt.savefig(f'{fname}.png', dpi=400)
plt.show()
my_nc.close()

# %%
pass
