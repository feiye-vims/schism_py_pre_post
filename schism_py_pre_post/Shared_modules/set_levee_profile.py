import shapefile
import numpy as np
from pylib import schism_grid, inside_polygon, proj
import pandas as pd
from hotstart_proc import nearest_neighbour
from download_nld import nld2map
import matplotlib.pyplot as plt


wdir = '/sciclone/scr10/feiye/NLD/v11.3/'
levee_name = 'LA_levees'

# read levee heights as xyz
_, xyz = nld2map(nld_fname=f'{wdir}/{levee_name}/System.geojson')
levee_x = xyz[:, 0]
levee_y = xyz[:, 1]
levee_height = xyz[:, 2]
levee_height[levee_height < 1] = 27
levee_height *= 0.3048
# plt.plot(np.sort(levee_height))
# plt.show()

proj(
    f'{wdir}/hgrid.utm.gr3', 0, 'epsg:26918',
    f'{wdir}/hgrid.ll', 0, 'epsg:4326',
)

gd = schism_grid(f'{wdir}/hgrid.utm.gr3')  # ; gd.save(f'{wdir}/hgrid.pkl')
gd_ll = schism_grid(f'{wdir}/hgrid.ll')  # ; gd.save(f'{wdir}/hgrid.pkl')
gd.lon = gd_ll.x
gd.lat = gd_ll.y


# find levee center line points in hgrid
shapefile_names = [f"{wdir}/la_levee_buffer.shp"]
ilevee = np.zeros(gd.dp.shape)
for shapefile_name in shapefile_names:
    sf = shapefile.Reader(shapefile_name)
    shapes = sf.shapes()
    for i, shp in enumerate(shapes):
        print(f'shp {i} of {len(shapes)}')
        poly_xy = np.array(shp.points).T
        ilevee += inside_polygon(np.c_[gd.x, gd.y], poly_xy[0], poly_xy[1])  # 1: true; 0: false
ilevee = ilevee.astype('bool')

gd.save(f'{wdir}/levees.gr3', value=ilevee)

II = nearest_neighbour(np.c_[gd.lon[ilevee], gd.lat[ilevee]], np.c_[levee_x, levee_y])
dist = np.sqrt((gd.lon[ilevee] - levee_x[II])**2 + (gd.lat[ilevee] - levee_y[II])**2)
short_dist = dist < 0.01  # degree, roughly 1000 m

# gd.dp[:] = 0
idx_levee_in_range = np.argwhere(ilevee)[:, 0][short_dist]
gd.dp[idx_levee_in_range] = - levee_height.astype(float)[II][short_dist]

gd.x = gd.lon
gd.y = gd.lat
gd.write_hgrid(f'{wdir}/levee_loaded.ll')

pass
