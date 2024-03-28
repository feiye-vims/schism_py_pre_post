import numpy as np
from scipy.spatial import KDTree
import geopandas as gpd
from pathlib import Path

from pylib_essentials.schism_file import read_schism_hgrid_cached
from pylib_essentials.utility_functions import inside_polygon
from matplotlib import pyplot as plt


# ---------------- input section ----------------
# chart data as a shapefile
# this is manually edited from the original navigation chart
# don't use the original chart data because the data points may be too sparse
# at places, e.g., in a channel with a maintained constant depth, there may be only
# one data point in the channel
sounding_shpfile = Path('/sciclone/schism10/Hgrid_projects/Charts/Savanna_Cooper/savannah_cooper_sounding.xyz_edited.shp')
# manually defined region in which to load the chart data
region_shpfile = Path('/sciclone/schism10/Hgrid_projects/Charts/Savanna_Cooper/manual_chart_region.shp')
crs_region = 'esri:102008'
# hgrid file
hg_file = Path('/sciclone/schism10/Hgrid_projects/Charts/Savanna_Cooper/hgrid.gr3')
# ----------------end input section ----------------

# Read in the data
sounding_data = gpd.read_file(sounding_shpfile)
# only keep the geometry column
sounding_data = sounding_data['geometry']
# extract xyz coordinates into an array
sounding_xyz = np.array([point.coords[0] for point in sounding_data])

# write the array to shapefile
# gpd.GeoDataFrame(geometry=gpd.points_from_xy(sounding_xyz[:,0], sounding_xyz[:,1], sounding_xyz[:,2])).to_file(
#    f"{sounding_shpfile.parent}/{sounding_shpfile.stem}_xyz{sounding_shpfile.suffix}")

# read the polygon shapefile
poly_shp = gpd.read_file(region_shpfile, crs=crs_region)
# shrink each polygon by 4 m
poly_shp['geometry'] = poly_shp['geometry'].buffer(-4)
# convert to lat/lon
poly_shp = poly_shp.to_crs('epsg:4326')

# extract the polygons to a list of nx2 numpy arrays
polygons = [np.array(poly.exterior.coords) for poly in poly_shp['geometry']]

hg = read_schism_hgrid_cached(hg_file)

# get the nodes in the polygons
in_channel = np.zeros_like(hg.dp, dtype=bool)
for polygon in polygons:
    idx = inside_polygon(np.c_[hg.x, hg.y], polygon[:, 0], polygon[:, 1]).astype(bool)
    in_channel[idx] = True

# for inpolygon nodes, find the nearest sounding_xyz point
# k-d tree
idx = KDTree(sounding_xyz[:, :2]).query(np.c_[hg.x, hg.y])[1]
nodal_depth_from_sounding_xyz = sounding_xyz[idx, 2]
dp_orig = hg.dp.copy()
hg.dp[in_channel] = nodal_depth_from_sounding_xyz[in_channel]
hg.dp = np.maximum(hg.dp, dp_orig)  # change the depth only if it is deeper than the original depth

# plot to check the changed nodes
plt.figure()
# add a 1:1 line for reference from the lower left to upper right
# Determine the range for the 1:1 line
min_val = min(min(dp_orig[in_channel]), min(hg.dp[in_channel]))
max_val = max(max(dp_orig[in_channel]), max(hg.dp[in_channel]))
# Add a 1:1 line for reference
plt.plot([min_val, max_val], [min_val, max_val], 'k--')
plt.scatter(dp_orig[in_channel], hg.dp[in_channel])
plt.axis('equal')
plt.show()

# output
hg.save(f"{hg_file.parent}/{hg_file.stem}_chart_loaded{hg_file.suffix}")

pass