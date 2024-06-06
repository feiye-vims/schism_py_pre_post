from pylib_essentials.schism_file import cread_schism_hgrid
import numpy as np

gd1 = cread_schism_hgrid('../v14a.gr3')
gd2 = cread_schism_hgrid('../hgrid.gr3')


diff = gd1.dp - gd2.dp
# find points with difference
idx = np.where(diff != 0)

# save difference to a file
np.savetxt('diff.txt', np.c_[gd1.x[idx], gd1.y[idx], diff[idx]], fmt='%f %f %f')
np.savetxt('gd2.txt', np.c_[gd2.x[idx], gd2.y[idx], gd2.dp[idx]], fmt='%f %f %f')
pass