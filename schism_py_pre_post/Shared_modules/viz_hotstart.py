from pylib import schism_grid
from hotstart_proc import Hotstart
import os
import numpy as np


my_hot = Hotstart(
    grid_info='/sciclone/schism10/feiye/STOFS3D-v4/Inputs/I21c/',
    hot_file='/sciclone/schism10/feiye/STOFS3D-v4/Inputs/I21c/hotstart.nc'
)

my_grid = schism_grid('/sciclone/schism10/feiye/STOFS3D-v4/Inputs/I21c/hgrid.gr3')

my_grid.plot(value=my_hot.tr_nd.val)
my_grid.dp[:] = my_hot.tr_nd.val[:, -1, 0]
my_grid.save('/sciclone/schism10/feiye/STOFS3D-v4/Inputs/I21c/surf_salt.gr3')

if os.path.exists(/sciclone/schism10/feiye/STOFS3D-v4/Inputs/I21c/hgrid.pkl):
    gd = schism_grid(/sciclone/schism10/feiye/STOFS3D-v4/Inputs/I21c/hgrid.pkl)
else:
    gd = schism_grid(/sciclone/schism10/feiye/STOFS3D-v4/Inputs/I21c/hgrid.gr3)
    gd.save(/sciclone/schism10/feiye/STOFS3D-v4/Inputs/I21c/hgrid.pkl)
pass