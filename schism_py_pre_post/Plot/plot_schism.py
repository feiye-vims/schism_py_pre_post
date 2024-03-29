from sympy import has_variety
import xarray
import numpy as np
import matplotlib.pyplot as plt
from pylib import schism_grid, grd2sms
import copy

'''

for i, file_id in enumerate(range(83, 91)):
    print(i)
    ds = xarray.open_dataset(f'/sciclone/scr10/feiye/STOFS3D-v5/RUN24a/outputs/out2d_{file_id}.nc')
    if i==0:
        elev = np.array(ds['elevation'])
    else:
        elev = np.concatenate((elev, np.array(ds['elevation'])), axis=0)
    ds.close()
#
# hg = schism_grid('/sciclone/schism10/feiye/STOFS3D-v5/Inputs/I24a/hgrid.pkl')  # hg.save('/sciclone/schism10/feiye/STOFS3D-v5/Inputs/I24a/hgrid.pkl')
# mask = elev[0, :] + hg.dp > 1e-6
# hg.dp = copy.deepcopy(elev[0, :])
# hg.dp[~mask] = -99999
# grd2sms(grd=hg, sms='/sciclone/schism10/feiye/STOFS3D-v5/Inputs/I24a/elev.2dm')

# ME_01049300 = (1778898, 1412496, 1397082, 1352813, 1367590)
# LA_08042522 = (578203,)
# LA_08030540 = (598585, 547975, 501620, 532565, 444206, 337719, 269403)
LA_MS = (267871, 296224, 349788)
for nd in LA_MS:
    plt.plot(elev[:, nd-1], label=f'{nd}')
    plt.legend()
plt.show()

pass
'''

import glob
output_dir = '/sciclone/home10/yinglong/DISKS/vims20/ChesapeakeBay/RUN112n_2010/outputs/'

# get dimensions
ds = xarray.open_dataset(f'{output_dir}/schout_1.nc')
nt, npoints = np.array(ds['WWM_1']).shape
rnday = 365
x = np.array(ds['SCHISM_hgrid_node_x'])
y = np.array(ds['SCHISM_hgrid_node_y'])

Hs_avg = np.zeros((npoints, ), dtype=float)
for i in range(rnday):
    output_file = f'{output_dir}/schout_{i+1}.nc'
    ds = xarray.open_dataset(output_file)
    print(f"Day {i+1} of {rnday}: avg Hs = {np.mean(np.array(ds['WWM_1'])[:])}")
    Hs_avg += np.mean(np.array(ds['WWM_1']), axis=0)/rnday
    ds.close()
np.savetxt('Hs_avg.txt', np.c_[x, y, Hs_avg])
pass