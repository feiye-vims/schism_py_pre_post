"""
Shift the datum of elev2D from xGEOID20b to NAVD88
"""


import numpy as np
from schism_py_pre_post.Utilities.util import vdatum_wrapper_pointwise, vdatum_preset
from pylib_experimental.schism_file import cread_schism_hgrid
import xarray as xr

wdir = '/sciclone/schism10/feiye/STOFS3D-v8/I14/InterpElev2D/'
gd = cread_schism_hgrid(f'{wdir}/fg.gr3')

x_bnd = gd.x[gd.iobn[0].astype(int)]
y_bnd = gd.y[gd.iobn[0].astype(int)]

z_shift = vdatum_wrapper_pointwise(
    x=x_bnd, y=y_bnd, z=np.zeros_like(x_bnd),
    conversion_para=vdatum_preset['xgeoid20b_to_navd88'],
    print_info='\nConverting from xgeoid20b to NAVD88:\n'
)

ds = xr.open_dataset(f'{wdir}/elev2D.th.nc')
ds['time_series'] = ds['time_series'] + z_shift[:, np.newaxis, np.newaxis]
ds.to_netcdf(f'{wdir}/elev2D.th.xgeoid20b_to_navd88.nc')
ds.close()
