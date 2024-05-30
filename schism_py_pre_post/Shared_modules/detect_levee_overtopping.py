import xarray as xr
import numpy as np
from pylib_essentials.schism_file import cread_schism_hgrid
from glob import glob

rundir = '/sciclone/schism10/feiye/STOFS3D-v7/Runs/R15e/'
levee_info_file = '/sciclone/schism10/feiye/STOFS3D-v7/Inputs/I14x/Bathy_edit/Levee/levee_top.gr3'

data = xr.open_mfdataset([
    f'{rundir}/outputs/out2d_29.nc',
    f'{rundir}/outputs/out2d_30.nc',
    f'{rundir}/outputs/out2d_31.nc',
    f'{rundir}/outputs/out2d_32.nc',
    f'{rundir}/outputs/out2d_33.nc',
    f'{rundir}/outputs/out2d_34.nc',
    f'{rundir}/outputs/out2d_37.nc',
])
elev = np.array(data['elevation'])
elev_max = np.max(elev, axis=0)

hgrid_levee = cread_schism_hgrid(levee_info_file)
ilevee = hgrid_levee.dp.copy().astype(bool)
levee_nd = np.argwhere(ilevee).flatten()

hgrid_obj = cread_schism_hgrid(f'{rundir}/hgrid.gr3')

dry = hgrid_obj.dp + elev_max < 1e-4
overtop_levee = ~dry & ilevee  # i.e, wet and on levee
if overtop_levee.any():  # dp is downward positive
    print(f'Levee overtopping detected on {sum(overtop_levee)} nodes.')

    with open(f'{rundir}/levee_overtop_nodes.txt', 'w') as f:
        f.write('lon lat z\n')
        for i in np.where(overtop_levee)[0]:
            f.write(f'{hgrid_obj.x[i]} {hgrid_obj.y[i]} {elev_max[i] + hgrid_obj.dp[i]}\n')
            print(f'Levee overtopped: node {i+1}, x={hgrid_obj.x[i]}, y={hgrid_obj.y[i]}, inundation depth={elev_max[i] + hgrid_obj.dp[i]}')

    np.savetxt(f'{rundir}/levee_overtop_nodes.txt', np.c_[hgrid_obj.x[overtop_levee], hgrid_obj.y[overtop_levee], elev_max[overtop_levee] + hgrid_obj.dp[overtop_levee]])

    # output elev_max with masked dry nodes
    elev_max[dry] = np.nan
    hgrid_obj.dp = elev_max
    hgrid_obj.grd2sms(f'{rundir}/maxelev.2dm')

pass
