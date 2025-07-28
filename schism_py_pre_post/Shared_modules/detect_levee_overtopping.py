"""
Detect levee overtopping based on the maximum water elevation on the levee top nodes
"""

import xarray as xr
import numpy as np
from schism_py_pre_post.Utilities.import_util import get_hgrid_reader
from glob import glob


# inputs
rundir = '/sciclone/schism10/feiye/STOFS3D-v8/R15d4_v7/'
# the levee top gr3 is written in Bathy_edit/Levee/set_levees
levee_top_gr3 = '/sciclone/schism10/feiye/STOFS3D-v8/I15g_v7/Bathy_edit/Levee/levee_top.gr3'
#  end inputs

# ------------------------------------read inputs------------------------------------
schism_grid = get_hgrid_reader()

hgrid_levee = schism_grid(levee_top_gr3)
ilevee = hgrid_levee.dp.copy().astype(bool)
levee_nd = np.argwhere(ilevee).flatten()

hgrid_obj = schism_grid(f'{rundir}/hgrid.gr3')
hgrid_obj.compute_ctr()

# read elevations
data = xr.open_mfdataset(glob(f'{rundir}/outputs/out2d_*.nc'))
elev = np.array(data['elevation'])
elev_max = np.max(elev, axis=0)

# ------------------------------------levee nodes-----------------------------------
# nodes on levee top
dry = hgrid_obj.dp + elev_max < 1e-1  # 0.1 m threshold
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
    
    print(f'maximum inundation depth on levee: {np.max(elev_max[overtop_levee] + hgrid_obj.dp[overtop_levee])}')

    # output elev_max with masked dry nodes
    elev_max[dry] = np.nan
    hgrid_obj.dp = elev_max
    hgrid_obj.grd2sms(f'{rundir}/maxelev.2dm')

# -----------------------------------levee elements-----------------------------------
# find elements on levee top (all element nodes on levee top)
ilevee = np.append(ilevee, [False, False])  # append two dummy nodes because elnode use -2 to indicate none-node
ilevee_ele = np.sum(ilevee[hgrid_obj.elnode], axis=1) >= 3
# remove the two dummy nodes
ilevee = ilevee[:-2]

# element is dry if any of its nodes is dry
dry = np.append(dry, [False, False])  # append two dummy nodes because elnode use -2 to indicate none-node
dry_ele = np.any(dry[hgrid_obj.elnode], axis=1)
dry = dry[:-2]  # remove the two dummy nodes

overtop_levee_ele = ilevee_ele & ~dry_ele
if overtop_levee_ele.any():
    print(f'Levee overtopping detected on {sum(overtop_levee_ele)} elements.')

    with open(f'{rundir}/levee_overtop_elements.txt', 'w', encoding='utf-8') as f:
        f.write('lon lat \n')
        for i in np.where(overtop_levee_ele)[0]:
            f.write(f'{hgrid_obj.xctr[i]} {hgrid_obj.yctr[i]} \n')
            print(f'Levee overtopped: element {i+1}, x={hgrid_obj.xctr[i]}, y={hgrid_obj.yctr[i]}')


print('Done.')
