"""
This script finds points within any regions (same format as those in Xmgredit)
"""


import numpy as np
from pylib import schism_grid, inside_polygon, read_schism_bpfile


def find_nodes_in_regions(xy, region_file_list):
    """
    Find points within any polygons specified in a shapefile.

    Inputs:
        xy: (np, 2) array of x and y coordinates.
        region_file_list: a list of strings (region file names).
    returns:
        i_inpoly: boolean array of size (np,) with True for points inside the regions.
        in_poly_idx: indices of the points inside the regions.
    """

    i_inpoly = np.zeros_like(xy[:, 0], dtype=bool)

    for region in region_file_list:
        bp = read_schism_bpfile(region, fmt=1)
        idx = inside_polygon(xy, bp.x, bp.y)
        i_inpoly = i_inpoly | idx

    in_poly_id = np.where(i_inpoly)[0]

    return i_inpoly, in_poly_id


def sample_usage():
    '''
    Sample usage of the function find_hgrid_nodes_regions
    '''
    # Load the hgrid
    gd = schism_grid('/sciclone/schism10/feiye/STOFS3D-v8/I15f_v7/Prop/hgrid.gr3')
    gd.compute_ctr()

    region_file_list = [
        '/sciclone/schism10/feiye/STOFS3D-v8/I15f_v7/Prop/nearshore_10m.rgn',
    ]

    # Find hgrid nodes within the regions
    i_inpoly, _ = find_nodes_in_regions(np.c_[gd.x, gd.y], region_file_list)

    print(f'Number of nodes inside the regions: {np.sum(i_inpoly)}')
    
    # do some operations with the nodes inside the polygons
    gd.save('/sciclone/schism10/feiye/STOFS3D-v8/I15f_v7/Prop/inear.gr3', value=i_inpoly.astype(int))


if __name__ == '__main__':
    sample_usage()
