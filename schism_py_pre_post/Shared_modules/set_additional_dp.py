from schism_py_pre_post.Geometry.inpoly import find_node_in_shpfiles
from pylib import schism_grid
import numpy as np

# v11.9, set upstream mississippi river levee height to 25 m, to prevent overtopping near the source
#        set Bonnet Carre Spill Way gate to 8.5 m
wdir = '/sciclone/schism10/feiye/STOFS3D-v4/Inputs/v11.9/'
gd_ll = schism_grid('/sciclone/schism10/feiye/STOFS3D-v4/Inputs/v11.7/v11.7.ll.pkl')  # gd_ll.save('/sciclone/schism10/feiye/STOFS3D-v4/Inputs/v11.7/v11.7.ll.pkl')
gd_utm = schism_grid('/sciclone/schism10/feiye/STOFS3D-v4/Inputs/v11.7/v11.7.utm.pkl')  # gd_utm.save('/sciclone/schism10/feiye/STOFS3D-v4/Inputs/v11.7/v11.7.utm.pkl')

for set_depth, shapefile in zip([-9, -20], [f'{wdir}/BonnetCarre_utm26918.shp', f'{wdir}/la_levee_center_line_upstream_missi_13m_buffer.shp']):
    i_inpoly = find_node_in_shpfiles(shapefile_names=[shapefile], gd=gd_utm)
    gd_ll.dp[i_inpoly] = np.minimum(gd_ll.dp[i_inpoly], set_depth)

gd_ll.save(f'{wdir}/hgrid.ll')