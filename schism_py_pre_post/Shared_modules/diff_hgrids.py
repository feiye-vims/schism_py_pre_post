"""
Check if two hgrids are different
"""

from copy import deepcopy
import numpy as np
from matplotlib import pyplot as plt
from pylib_experimental.schism_file import cread_schism_hgrid

hgrid1 = '/sciclone/schism10/feiye/STOFS3D-v8/I15a_v7/Bathy_edit/hgrid_dem_edit.ll'
hgrid2 = '/sciclone/schism10/feiye/STOFS3D-v8/I15a_v7/Bathy_edit2/hgrid_dem_edit.ll'

gd1 = cread_schism_hgrid(hgrid1)
gd2 = cread_schism_hgrid(hgrid2)

assert gd1.np == gd2.np, 'Number of nodes are different'
assert gd1.ne == gd2.ne, 'Number of elements are different'
assert np.allclose(gd1.x, gd2.x), 'x coordinates are different'
assert np.allclose(gd1.y, gd2.y), 'y coordinates are different'

dp_diff = np.abs(gd1.dp - gd2.dp)
hgrid_diff = deepcopy(gd1)
hgrid_diff.dp = gd2.dp - gd1.dp
hgrid_diff.plot(fmt=1, clim=[-1, 1])
plt.show()
print(
    'Maximum difference in depths:\n'
    f'{np.max(abs(dp_diff))} meters\n'
    f'{np.max(abs(dp_diff)) * 100} cm\n'
    f'{np.max(abs(dp_diff)) * 1000} mm\n'
    '\n'
    'Mean difference in depths:\n'
    f'{np.mean(abs(dp_diff))} meters\n'
    f'{np.mean(abs(dp_diff)) * 100} cm\n'
    f'{np.mean(abs(dp_diff)) * 1000} mm\n'
)
assert np.allclose(gd1.dp, gd2.dp), 'depths are different'
print('Two hgrids are the same')