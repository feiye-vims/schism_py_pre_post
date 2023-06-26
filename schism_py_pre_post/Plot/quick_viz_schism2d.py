try:
    from schism_py_pre_post.Grid.Hgrid_extended import read_schism_hgrid_cached, read_schism_vgrid_cached
except ImportError:
    from spp_essentials.Hgrid_extended import read_schism_hgrid_cached, read_schism_vgrid_cached
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from IPython.display import set_matplotlib_formats
set_matplotlib_formats('svg')


gd_fname = '/sciclone/schism10/feiye/STOFS3D-v6/Inputs/I15/hgrid.gr3'
gd = read_schism_hgrid_cached(gd_fname, overwrite_cache=False)

vg_fname = '/sciclone/schism10/feiye/STOFS3D-v6/Inputs/I15/vgrid.in'
vg = read_schism_vgrid_cached(vg_fname, overwrite_cache=False)

var_name = "salinity"
fname = f'/sciclone/home/feiye/Sync/{var_name}_321.nc'
# fname = '/sciclone/scr10/lcui01/ICOGS3D/outputs_RUN13u6/temperature_90.nc'

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

pass
