from pylib_essentials.schism_file import grd2sms, read_schism_hgrid_cached, read_schism_vgrid_cached
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from IPython.display import set_matplotlib_formats
from pathlib import Path
import xarray as xr
import copy
set_matplotlib_formats('svg')


gd_fname = '/sciclone/schism10/feiye/Test_Project/Runs/R98a/hgrid.gr3'
gd = read_schism_hgrid_cached(gd_fname, overwrite_cache=False)

vg_fname = '/sciclone/schism10/feiye/Test_Project/Runs/R98a/vgrid.in'
vg = read_schism_vgrid_cached(vg_fname, overwrite_cache=False)

var_name = "elevation"
fnames = [f'/sciclone/schism10/feiye/Test_Project/Runs/R99a/outputs/out2d_1.nc']
my_nc = xr.open_mfdataset(fnames)

it = -1
isurf = True
caxis = [0, 100]
xlim = None  # [-77, -75]
ylim = None  # [37, 40]
plt.figure(figsize=(7, 7))

var = np.array(my_nc[var_name])

if isurf:
    if len(var.shape) == 3:
        value = var[it, :, -1]
    elif len(var.shape) == 2:
        value = np.max(var[:, :], axis=0)
else:
    value = var[it, :, :]
    value = value[np.arange(gd.np), vg.kbp]

disturbance = copy.deepcopy(value)
land = gd.dp < 0.0
disturbance[land] = value[land] + gd.dp[land]
idx = disturbance > 10
np.savetxt(f'{Path(fnames[0]).parent}/high_disturbance.xyz', np.vstack((gd.x[idx], gd.y[idx], disturbance[idx])).T, fmt='%.6f', delimiter=' ')

gd.dp = disturbance
grd2sms(gd, str(Path(fnames[0]).with_suffix('.disturbance.2dm')))

gd.dp = value
grd2sms(gd, str(Path(fnames[0]).with_suffix(f'.{var_name}.2dm')))

gd.plot_grid(fmt=1, value=value, clim=caxis, levels=31, ticks=np.arange(caxis[0], caxis[1], 2), cmap='jet', xlim=xlim, ylim=ylim)
plt.gca().set_aspect('equal', 'box')
plt.savefig(f'{fnames[0]}.png', dpi=400)
plt.show()
my_nc.close()

pass
