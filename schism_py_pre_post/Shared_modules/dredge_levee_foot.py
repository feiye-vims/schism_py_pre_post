# %%
from schism_py_pre_post.Geometry.inpoly import find_node_in_shpfiles
from pylib import schism_grid, grd2sms
import numpy as np


def row_unique(a):
    unique = np.sort(a)
    duplicates = unique[:,  1:] == unique[:, :-1]
    unique[:, 1:][duplicates] = -2
    return unique

# %%
gd = schism_grid('/sciclone/schism10/feiye/STOFS3D-v4/Inputs/I23o/Dredge_levee_foot/hgrid.utm.gr3')

# %%
i_levee_foot = find_node_in_shpfiles(gd=gd, shapefile_names=['/sciclone/schism10/feiye/STOFS3D-v4/Inputs/I23o/Dredge_levee_foot/Shapefiles/levee_foot.shp']).astype(bool)

# gd.dp[:] = 0
# gd.dp[i_levee_foot] = 1
# grd2sms(gd, '/sciclone/schism10/feiye/STOFS3D-v4/Inputs/I23o/Dredge_levee_foot/levee_foot.2dm')

# find neighboring nodes of each foot
gd.compute_all()

ine_levee_foot = gd.ine[i_levee_foot, :]
n_levee_foot, max_nei = ine_levee_foot.shape
inp_levee_foot = row_unique(np.reshape(gd.elnode[ine_levee_foot, :], (n_levee_foot, -1)))
dp_padded = np.r_[gd.dp, np.array([-9999, -9999])]  # pad a high elevation for index "-2" in inp
dredged_levee_foot_dp = np.max(dp_padded[inp_levee_foot], axis=1)

gd.dp[i_levee_foot] = dredged_levee_foot_dp[:]

# output
gd.save('/sciclone/schism10/feiye/STOFS3D-v4/Inputs/I23o/Dredge_levee_foot/dredged.utm.gr3')

gd.x, gd.y = gd.proj(prj0='epsg:26918', prj1='epsg:4326')
gd.save('/sciclone/schism10/feiye/STOFS3D-v4/Inputs/I23o/Dredge_levee_foot/dredged.ll')

grd2sms(gd, '/sciclone/schism10/feiye/STOFS3D-v4/Inputs/I23o/Dredge_levee_foot/dredged.2dm')