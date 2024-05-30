#!/usr/bin/env python3

'''
This script converts a shapefile to a region file for use in SCHISM
'''

from pathlib import Path

import numpy as np
import geopandas as gpd

from pylib import schism_bpfile
from pylib_experimental.schism_file import cread_schism_hgrid

wdir = '/sciclone/schism10/feiye/STOFS3D-v7/Inputs/I18c/Prop/'
out_name = 'iest'
# input_filename = '/sciclone/schism10/feiye/STOFS3D-v7/Inputs/I18c/Prop/nearshore_10m.rgn'
input_filename = '/sciclone/schism10/feiye/STOFS3D-v7/Inputs/I18c/Prop/R18c_tvd_polygons.shp'
if Path(input_filename).suffix in ['.rgn', 'reg']:
    # convert a region file to a shapefile
    region = schism_bpfile()
    region.read(fname=input_filename)
    from shapely.geometry import Polygon
    poly = Polygon(np.c_[region.x, region.y])
    polygons_gdf = gpd.GeoDataFrame(geometry=[poly], crs='epsg:4326')
elif Path(input_filename).suffix in ['.shp']:
    polygons_gdf = gpd.read_file(input_filename)
else:
    raise ValueError('Input file must be either a region file or a shapefile')

# buffer each polygon by 1e-4 degrees (about 10 meters)
crs = polygons_gdf.crs.to_string()
polygons_gdf = gpd.GeoDataFrame(geometry=polygons_gdf.to_crs('esri:102008').buffer(10).to_crs(crs))

gd = cread_schism_hgrid('/sciclone/schism10/feiye/STOFS3D-v7/Inputs/I18c/hgrid.gr3')
# make a geopandas dataframe with the same crs as the shapefile
gd_points_gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy(gd.x, gd.y), crs='epsg:4326')

# intersect gd_points and shp to find points inside polygons
in_poly = gpd.sjoin(gd_points_gdf, polygons_gdf, how="inner", predicate='within').index
# write to *.gr3
idx = np.zeros(gd.np, dtype=int)
idx[in_poly] = 1
gd.save(f'{wdir}/{out_name}.gr3', value=idx)

pass