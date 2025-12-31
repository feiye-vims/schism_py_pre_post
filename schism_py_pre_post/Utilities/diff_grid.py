from pylib_experimental.schism_file import cread_schism_hgrid as read_schism_hgrid
from copy import deepcopy
import numpy as np
from pylib import grd2sms  #, read_schism_hgrid

gd1 = read_schism_hgrid('/sciclone/schism10/feiye/STOFS3D-v8/I15y_v7/Set_Levee/hgrid_with_levees_xGEOID.gr3')
gd2 = read_schism_hgrid('/sciclone/schism10/feiye/STOFS3D-v8/I15w_v7/hgrid_levee_2025_xgeoid20b.gr3')

is_equal = np.allclose(gd1.dp, gd2.dp, atol=1e-6, rtol=1e-6)

# find points with difference in x, y, or dp
if is_equal:
    print("The grids are equal.")
else:
    print("The grids are not equal.")
    # diffx = gd1.x - gd2.x
    # diffy = gd1.y - gd2.y
    diff_dp = gd1.dp - gd2.dp
    diff_mask = np.abs(diff_dp) > 1e-4  # 0.1 mm tolerance

idx = np.where(diff_mask)[0]

# save the diff as a new grid dp
gd_diff = deepcopy(gd1)
gd_diff.dp = diff_mask.astype(int)

# output
out_dir = '/sciclone/schism10/feiye/STOFS3D-v8/I15y_v7/'
grd2sms(gd_diff, f'{out_dir}/diff.2dm')

# save difference to a file
np.savetxt(f'{out_dir}/diff_gd1.txt', np.c_[gd1.x[idx], gd1.y[idx], diff_dp[idx]], fmt='%f %f %f')
np.savetxt(f'{out_dir}/diff_gd2.txt', np.c_[gd2.x[idx], gd2.y[idx], gd2.dp[idx]], fmt='%f %f %f')
pass