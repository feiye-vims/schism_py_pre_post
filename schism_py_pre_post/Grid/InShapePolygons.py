"""
This script finds points within any polygons specified in a shapefile.
"""


import numpy as np
import geopandas as gpd
from pylib import schism_grid


def find_nodes_in_polygons(xy, shapefile):
    """
    Find points within any polygons specified in a shapefile.

    Inputs:
        xy: (np, 2) array of x and y coordinates.
        shapefile: a string (shapefile name) or a GeoDataFrame.
    returns:
        i_inpoly: boolean array of size (np,) with True for points inside the polygons.
        in_poly_idx: indices of the points inside the polygons.
    """

    if shapefile is str:
        gdf = gpd.read_file(shapefile)
    elif isinstance(shapefile, gpd.GeoDataFrame):
        gdf = shapefile
    else:
        raise ValueError('shapefile must be a string or a GeoDataFrame')

    xy_gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy(xy[:, 0], xy[:, 1]))
    joined_gdf = gpd.sjoin(xy_gdf, gdf, how="inner", predicate='within')

    # Get the indices of the points inside the polygons
    in_poly_id = joined_gdf.index.to_numpy()
    i_inpoly = np.zeros_like(xy[:, 0], dtype=bool)
    i_inpoly[in_poly_id] = True

    return i_inpoly, in_poly_id


def sample_usage():
    '''
    Sample usage of the function find_hgrid_nodes_in_polygons
    '''
    # Load the hgrid
    gd = schism_grid('/sciclone/schism10/feiye/STOFS3D-v8/I15f_v7/Prop/hgrid.gr3')
    gd.compute_ctr()

    # Load the shapefile
    gdf = gpd.read_file('/sciclone/schism10/feiye/STOFS3D-v8/I15f_v7/Prop/coastal_tvd.shp')
    gdf = gdf.to_crs('epsg:4326')

    # Find hgrid nodes within the polygons
    i_inpoly, _ = find_nodes_in_polygons(np.c_[gd.xctr, gd.yctr], gdf)

    print(f'Number of nodes inside the polygons: {np.sum(i_inpoly)}')
    
    # do some operations with the nodes inside the polygons
    gd.write_prop("in_poly.prop", i_inpoly.astype(int), fmt='{:d}')
    upwind_ele_idx = ~i_inpoly
    upwind_nodes = np.unique(gd.elnode[upwind_ele_idx].reshape(-1, ))
    upwind_nodes = upwind_nodes[upwind_nodes >= 0]  # remove dummy nodes of value -2
    tvd_nodes = np.ones((gd.np, ), dtype=int)
    tvd_nodes[upwind_nodes] = 0
    gd.save('/sciclone/schism10/feiye/STOFS3D-v8/I15f_v7/Prop/iest.gr3', value=tvd_nodes)


if __name__ == '__main__':
    sample_usage()
