import numpy as np
from pylib import inside_polygon
from scipy.spatial import cKDTree
import shapefile


def find_ele_node_in_shpfile(shapefile_name, grid):
    '''
    Find element/node index within one or more polygons defined in a shapefile
        shapefile_name: file contains polygon(s)
        grid: schism_grid instance
    '''
    grid.compute_ctr()

    sf = shapefile.Reader(shapefile_name)
    shapes = sf.shapes()

    ele_ind_list = []
    for shp in shapes:
        poly_xy = np.array(shp.points).T
        ind = inside_polygon(np.c_[grid.xctr, grid.yctr], poly_xy[0], poly_xy[1])  # 1: true; 0: false
        ind = ind.astype('bool')
        ele_ind_list.append(ind)

    node_ind_list = []
    for shp in shapes:
        poly_xy = np.array(shp.points).T
        ind = inside_polygon(np.c_[grid.x, grid.y], poly_xy[0], poly_xy[1])  # 1: true; 0: false
        ind = ind.astype('bool')
        node_ind_list.append(ind)

    return [ele_ind_list, node_ind_list]
