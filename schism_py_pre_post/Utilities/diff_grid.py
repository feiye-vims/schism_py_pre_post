import os
from pylib_experimental.schism_file import cread_schism_hgrid as read_schism_hgrid
from copy import deepcopy
import numpy as np
from pylib import grd2sms  #, read_schism_hgrid

gd1 = read_schism_hgrid('/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v32/v32b.gr3')
gd2 = read_schism_hgrid('/sciclone/schism10/hjyoo/task/task10_Atlantic/RUN100b/hgrid.gr3')
output_dir = '/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v32/diff'

os.makedirs(output_dir, exist_ok=True)

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

# historgram of differences
print(f"Number of points with dp difference > 0.1 mm: {np.sum(diff_mask)}")
print(f"Max dp difference: {np.max(np.abs(diff_dp))} m")

n_largest = 100
largest_diff_idx = np.argsort(np.abs(diff_dp))[::-1][:n_largest]
print(f"Top {n_largest} largest dp differences:")
for idx in largest_diff_idx:
    print(f"Node {idx+1}: lat, lon = ({gd1.y[idx]}, {gd1.x[idx]}), dp1 = {gd1.dp[idx]}, dp2 = {gd2.dp[idx]}, diff = {diff_dp[idx]} m")
np.savetxt(f'{output_dir}/largest_diff_dp.txt', np.c_[gd1.x[largest_diff_idx], gd1.y[largest_diff_idx], gd1.dp[largest_diff_idx], gd2.dp[largest_diff_idx], diff_dp[largest_diff_idx]], fmt='%f %f %f %f %f', header='lon lat dp1 dp2 diff_dp')

import matplotlib.pyplot as plt
plt.hist(diff_dp, bins=50)
plt.xlabel('Difference in dp (m)')
plt.ylabel('Frequency')
plt.title('Histogram of dp differences between two grids')
plt.show()

idx = np.where(diff_mask)[0]

# save the diff as a new grid dp
gd_diff = deepcopy(gd1)
gd_diff.dp = diff_mask.astype(int)
grd2sms(gd_diff, f'{output_dir}/diff_mask.2dm')
gd_diff.dp = diff_dp
grd2sms(gd_diff, f'{output_dir}/diff_dp.2dm')

# save difference to a file
np.savetxt(f'{output_dir}/diff_gd1.txt', np.c_[gd1.x[idx], gd1.y[idx], diff_dp[idx]], fmt='%f %f %f')
np.savetxt(f'{output_dir}/diff_gd2.txt', np.c_[gd2.x[idx], gd2.y[idx], gd2.dp[idx]], fmt='%f %f %f')
pass