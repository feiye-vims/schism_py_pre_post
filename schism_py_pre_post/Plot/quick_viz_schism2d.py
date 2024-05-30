from pylib_experimental.schism_file import cread_schism_hgrid, read_schism_vgrid_cached
from pylib import grd2sms, schism_grid
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from IPython.display import set_matplotlib_formats
from pathlib import Path
import xarray as xr
import copy
set_matplotlib_formats('svg')

print('reading data ...\n')

gd_fname = '/sciclone/schism10/feiye/STOFS3D-v6/Inputs/v6.1.static_inputs_snapshot20231105/hgrid.gr3'
gd = cread_schism_hgrid(gd_fname)

vg_fname = '/sciclone/schism10/feiye/STOFS3D-v6/Runs/RUN24a/vgrid.in'
vg = read_schism_vgrid_cached(vg_fname, overwrite_cache=False)

var_name = "salt_surface"
fnames = [f'/sciclone/home/feiye/TEMP/v2.1/stofs_3d_atl.t12z.field2d_f025_036.nc']
my_nc = xr.open_mfdataset(fnames)

elements = my_nc['SCHISM_hgrid_face_nodes'].values
elements.shape
i3_idx = np.argwhere(gd.i34==3).flatten()
np.array_equal(gd.elnode[i3_idx, :3], elements[i3_idx, :3]-1)

it = -1
isurf = True
caxis = [0, 30]
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

# output high disturbance points
# idx = disturbance > 10
# np.savetxt(f'{Path(fnames[0]).parent}/high_disturbance.xyz', np.vstack((gd.x[idx], gd.y[idx], disturbance[idx])).T, fmt='%.6f', delimiter=' ')

i_positive_disturbance = disturbance > 0  # i.e., abnormal inundation
disturbance[~i_positive_disturbance] = -9999  # set negative disturbance to -9999
gd.dp = disturbance
grd2sms(gd, str(Path(fnames[0]).with_suffix('.positive_disturbance.2dm')))

gd.dp = value
grd2sms(gd, str(Path(fnames[0]).with_suffix(f'.{var_name}.2dm')))

gd.plot_grid(fmt=1, value=value, clim=caxis, levels=31, ticks=np.arange(caxis[0], caxis[1], 2), cmap='jet', xlim=xlim, ylim=ylim)
plt.gca().set_aspect('equal', 'box')
plt.savefig(f'{fnames[0]}.png', dpi=400)
plt.show()
my_nc.close()

pass
