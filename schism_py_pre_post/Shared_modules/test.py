from pylib import harmonic_analysis, schism_grid
import numpy as np
from pylib_essentials.schism_file import read_schism_hgrid_cached
from shapely.geometry import Polygon
import geopandas as gpd
from pathlib import Path

filename = Path('/sciclone/schism10/Hgrid_projects/DEMs/hgrid.dem_id.2dm')
hg = read_schism_hgrid_cached(filename)
nodes = np.c_[hg.x, hg.y, hg.dp]
nodes = np.r_[nodes, np.c_[np.nan, np.nan, np.nan]]  # add a dummy node to accomodate for -1 in elnode

# replace -2 with -1
hg.elnode[hg.elnode == -2] = -1

elnode_dp = nodes[:, 2][hg.elnode]
for n in [3, 4]:  # triangle and quad
    idx = np.argwhere(hg.i34 == n).flatten()
    valid = np.all(elnode_dp[idx, 1:n] == elnode_dp[idx, :n-1], axis=1)  # check if all nodes have the same z
    elnode_dp[idx[~valid], :] = -999  # set invalid elements to -999, because the DEM interpolation is node based
dpe = elnode_dp[:, 0].astype(int)

polygons = []
for i, [i34, elnode] in enumerate(zip(hg.i34, hg.elnode)):
    elnode = elnode[:i34]
    element_coords = nodes[elnode, :2]
    polygons.append(Polygon(element_coords))

gdf = gpd.GeoDataFrame({'geometry': polygons, 'z': dpe})
gdf = gdf.dissolve(by='z')
gdf.to_file(filename.with_suffix('.shp'))

# Convert coordinates to Polygons
# %%
import pygrib
import copy
import numpy as np

filename = '/sciclone/schism10/feiye/STOFS3D-v4/AWIPS/estofs.t00z.conus.east.cwl.grib2'
outfile = '/sciclone/schism10/feiye/STOFS3D-v4/AWIPS/test.grib2'

gr = pygrib.open(filename)
with open(outfile, 'wb') as outgrb:
    for g in gr:
        # print(g.typeOfLevel, g.level, g.name, g.validDate, g.analDate, g.forecastTime)
        # vals = copy.deepcopy(g['values'])
        # vals[vals.mask == False] = 0
        g['values'] -= g['values']
        g['values'] -= g['values']
        g['values'] = g['values'] + 10.0
        msg = g.tostring()
        outgrb.write(msg)
outgrb.close()
gr.close()
pass

# %%
import xarray as xr
from cfgrib.xarray_to_grib import to_grib


filename = '/sciclone/schism10/feiye/STOFS3D-v4/AWIPS/estofs.t00z.conus.east.cwl.grib2'
ds=xr.open_dataset(filename, engine='cfgrib')

ds.to_netcdf('test.nc')
 

# %%
from pyproj import Proj, transform
inProj = Proj('esri:102009')
outProj = Proj('epsg:4326')

x2,y2 = transform(inProj,outProj, 2681.9,-263.8)

pass
