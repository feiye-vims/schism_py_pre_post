"""
Conditionally transfer hgrid's xy from one to another
"""

from copy import deepcopy
import numpy as np
from matplotlib import pyplot as plt
from pylib_experimental.schism_file import cread_schism_hgrid

hgrid1 = '/sciclone/schism10/feiye/STOFS3D-v8/R15e_v7/hgrid.gr3'
hgrid2 = '/sciclone/schism10/feiye/STOFS3D-v8/R15e_v7/hgrid_edit.gr3'

gd1 = cread_schism_hgrid(hgrid1)
gd2 = cread_schism_hgrid(hgrid2)

assert gd1.np == gd2.np, 'Number of nodes are different'
assert gd1.ne == gd2.ne, 'Number of elements are different'

nodal_dist = np.sqrt((gd1.x - gd2.x) ** 2 + (gd1.y - gd2.y) ** 2)
significant_diff = nodal_dist > 1e-5
print(f'Number of nodes with significant difference: {np.sum(significant_diff)}')

gd1.x[significant_diff] = gd2.x[significant_diff]
gd1.y[significant_diff] = gd2.y[significant_diff]

gd1.save('/sciclone/schism10/feiye/STOFS3D-v8/R15e_v7/hgrid_xy_transferred.ll')
print('Done.')