from pylib_experimental.schism_file import cread_schism_hgrid
from pylib import grd2sms
import numpy as np

gd1 = cread_schism_hgrid('/sciclone/schism10/feiye/STOFS3D-v7/Inputs/v7_static_inputs_20240911/hgrid_close_MTG_IWW.gr3')
gd2 = cread_schism_hgrid('/sciclone/schism10/feiye/STOFS3D-v8/I13a/Bathy_edit/hgrid_dem_edit.ll')

is_equal = np.allclose(gd1.dp, gd2.dp, atol=1e-6, rtol=1e-6)

diff = gd1.dp - gd2.dp
# find points with difference
idx = np.where(diff != 0)

# save the diff as a new grid
gd1.dp = diff
grd2sms(gd1, '/sciclone/schism10/feiye/STOFS3D-v8/I13a/Bathy_edit/v7-13a.2dm')

# save difference to a file
np.savetxt('diff.txt', np.c_[gd1.x[idx], gd1.y[idx], diff[idx]], fmt='%f %f %f')
np.savetxt('gd2.txt', np.c_[gd2.x[idx], gd2.y[idx], gd2.dp[idx]], fmt='%f %f %f')
pass