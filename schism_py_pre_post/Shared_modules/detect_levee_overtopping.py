import xarray as xr
import numpy as np
from pylib_essentials.schism_file import cread_schism_hgrid
from pylib import schism_grid
from glob import glob

rundir = '/sciclone/schism10/feiye/STOFS3D-v7/Runs/Rv7_Francine_update/'
levee_info_file = '/sciclone/schism10/feiye/STOFS3D-v7/Inputs/I14x/Bathy_edit/Levee/levee_top.gr3'

data = xr.open_mfdataset([
    f'{rundir}/outputs/out2d_1.nc',
    f'{rundir}/outputs/out2d_2.nc',
    f'{rundir}/outputs/out2d_3.nc',
])
elev = np.array(data['elevation'])
elev_max = np.max(elev, axis=0)

hgrid_levee = schism_grid(levee_info_file)
ilevee = hgrid_levee.dp.copy().astype(bool)
levee_nd = np.argwhere(ilevee).flatten()

hgrid_obj = schism_grid(f'{rundir}/hgrid.gr3')

dry = hgrid_obj.dp + elev_max < 1e-1
overtop_levee = ~dry & ilevee  # i.e, wet and on levee
if overtop_levee.any():  # dp is downward positive
    print(f'Levee overtopping detected on {sum(overtop_levee)} nodes.')

    with open(f'{rundir}/levee_overtop_nodes.txt', 'w', encoding='utf-8') as f:
        f.write('lon lat z\n')
        for i in np.where(overtop_levee)[0]:
            f.write(f'{hgrid_obj.x[i]} {hgrid_obj.y[i]} {elev_max[i] + hgrid_obj.dp[i]}\n')
            print(
                f'Levee overtopped: node {i+1}, x={hgrid_obj.x[i]}, y={hgrid_obj.y[i]},'
                f'inundation depth={elev_max[i] + hgrid_obj.dp[i]}'
            )

    # output elev_max with masked dry nodes
    elev_max[dry] = np.nan
    hgrid_obj.dp = elev_max
    hgrid_obj.grd2sms(f'{rundir}/maxelev.2dm')

pass
