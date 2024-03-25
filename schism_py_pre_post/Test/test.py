import geopandas as gpd
import numpy as np
from pathlib import Path

from scipy.spatial import KDTree
from matplotlib import pyplot as plt

from pylib_essentials.schism_file import read_schism_hgrid_cached, grd2sms
from pylib_essentials.utility_functions import inside_polygon

hg_file = Path('/sciclone/schism10/Hgrid_projects/Charts/Mississippi/hgrid.gr3')

# prepare the sounding data
sounding_shpfile = Path('/sciclone/schism10/Hgrid_projects/Charts/Mississippi/chart_merged_points.shp')

sounding_data = gpd.read_file(sounding_shpfile)
# only keep the geometry column
sounding_data = sounding_data['geometry']
# extract xyz coordinates into an array
sounding_xyz = np.array([point.coords[0] for point in sounding_data])

# create a new geodataframe with the xyz coordinates, and set z as a column
sounding_gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy(sounding_xyz[:,0], sounding_xyz[:,1]))
sounding_gdf['z'] = sounding_xyz[:,2]

# write gdf to shapefile
sounding_gdf.to_file(f"{sounding_shpfile.parent}/{sounding_shpfile.stem}_xyz{sounding_shpfile.suffix}")

# prepare the polygon data
channel_shpfile = Path('/sciclone/schism10/Hgrid_projects/Charts/Mississippi/chart_region_from_sms_quads.shp', crs='esri:102008')
channel_polys = gpd.read_file(channel_shpfile)

# enlarge the polygons by 4 m
channel_polys['geometry'] = channel_polys['geometry'].buffer(4)

# convert to lat/lon
channel_polys = channel_polys.to_crs('epsg:4326')

# extract the polygons to a list of nx2 numpy arrays
polygons = [np.array(poly.exterior.coords) for poly in channel_polys['geometry']]


hg = read_schism_hgrid_cached(hg_file)

# get the nodes in the polygons
in_channel = np.zeros_like(hg.dp, dtype=bool)
for polygon in polygons:
    idx = inside_polygon(np.c_[hg.x, hg.y], polygon[:, 0], polygon[:, 1]).astype(bool)
    in_channel[idx] = True
# diagnostic plot
plt.figure()
hg.plot(value=in_channel.astype(int), fmt=1)

# for inpolygon nodes, find the nearest sounding_xyz point
# k-d tree
idx = KDTree(sounding_xyz[:, :2]).query(np.c_[hg.x, hg.y])[1]
dp_sounding = sounding_xyz[idx, 2]


# get the nodes in the NCF (National Channel Framework) polygons
NCF_shpfile = Path('/sciclone/schism10/Hgrid_projects/Charts/Mississippi/channel_quarter_NCF.shp')
NCF_data = gpd.read_file(NCF_shpfile)

# remove the polygons that are not in the bounding box of the hgrid
NCF_data = NCF_data.cx[hg.x.min():hg.x.max(), hg.y.min():hg.y.max()]
# convert any multipolygons to polygons
NCF_data = NCF_data.explode()
# enlarge the polygons by 4 m
NCF_data['geometry'] = NCF_data['geometry'].to_crs('esri:102008').buffer(4).to_crs('epsg:4326')
# extract the polygons to a list of nx2 numpy arrays
NCF_polygons = [np.array(poly.exterior.coords) for poly in NCF_data['geometry']]

# get the nodes in the polygons
dp_NCF = np.zeros_like(hg.dp, dtype=float)
in_NCF = np.zeros_like(hg.dp, dtype=bool)
for maintained_depth, polygon in zip(NCF_data['depthmaint'], NCF_polygons):
    idx = inside_polygon(np.c_[hg.x, hg.y], polygon[:, 0], polygon[:, 1]).astype(bool)
    in_NCF[idx] = True
    dp_NCF[idx] = maintained_depth * 0.3048  # convert from feet to meters
# diagnostic plot
plt.figure()
hg.plot(value=in_NCF.astype(int), fmt=1)

dp_orig = hg.dp.copy()
hg.dp = np.maximum(dp_orig, dp_NCF, dp_sounding)  # change the depth only if it is deeper than the original depth

# plot to check the changed nodes
plt.figure()
idx = abs(dp_orig - hg.dp) > 1e-5
# add a 1:1 line for reference from the lower left to upper right
# Determine the range for the 1:1 line
min_val = min(min(dp_orig[idx]), min(hg.dp[idx]))
max_val = max(max(dp_orig[idx]), max(hg.dp[idx]))
# Add a 1:1 line for reference
plt.plot([min_val, max_val], [min_val, max_val], 'k--')
plt.scatter(dp_orig[idx], hg.dp[idx])
plt.axis('equal')
plt.show()

# output
grd2sms(hg, (f"{hg_file.parent}/{hg_file.stem}_chart_loaded.2dm"))

pass
