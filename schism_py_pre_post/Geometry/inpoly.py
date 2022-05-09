import numpy as np
from pylib import inside_polygon
from scipy.spatial import cKDTree
import shapefile


def find_pts_in_shpfiles(shapefile_names=[], pts=None):
    '''
    Find node index within one or more polygons defined in a shapefile
        shapefile_name: file contains polygon(s)
        gd: schism_grid instance
    '''
    pts = np.array(pts)
    idx = np.zeros(pts.shape[0])
    for shapefile_name in shapefile_names:
        sf = shapefile.Reader(shapefile_name)
        shapes = sf.shapes()
        for i, shp in enumerate(shapes):
            print(f'shp {i+1} of {len(shapes)}')
            poly_xy = np.array(shp.points).T
            idx += inside_polygon(np.c_[pts[:, 0], pts[:, 1]], poly_xy[0], poly_xy[1])  # 1: true; 0: false
    idx = idx.astype('bool')
    return idx

def find_ele_in_shpfiles(shapefile_names=[], gd=None):
    '''
    Find node index within one or more polygons defined in a shapefile
        shapefile_name: file contains polygon(s)
        gd: schism_grid instance
    '''
    gd.compute_ctr()
    idx = np.zeros(gd.dpe.shape)
    for shapefile_name in shapefile_names:
        sf = shapefile.Reader(shapefile_name)
        shapes = sf.shapes()
        for i, shp in enumerate(shapes):
            print(f'shp {i+1} of {len(shapes)}')
            poly_xy = np.array(shp.points).T
            idx += inside_polygon(np.c_[gd.xctr, gd.yctr], poly_xy[0], poly_xy[1])  # 1: true; 0: false
    idx = idx.astype('bool')
    return idx


def find_node_in_shpfiles(shapefile_names=[], gd=None):
    '''
    Find element index within one or more polygons defined in a shapefile
        shapefile_name: file contains polygon(s)
        gd: schism_grid instance
    '''
    idx = np.zeros(gd.dp.shape)
    for shapefile_name in shapefile_names:
        sf = shapefile.Reader(shapefile_name)
        shapes = sf.shapes()
        for i, shp in enumerate(shapes):
            print(f'shp {i+1} of {len(shapes)}')
            poly_xy = np.array(shp.points).T
            idx += inside_polygon(np.c_[gd.x, gd.y], poly_xy[0], poly_xy[1])  # 1: true; 0: false
    idx = idx.astype('bool')
    return idx


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
    for i, shp in enumerate(shapes):
        print(f'shp {i+1} of {len(shapes)}')
        poly_xy = np.array(shp.points).T
        ind = inside_polygon(np.c_[grid.xctr, grid.yctr], poly_xy[0], poly_xy[1])  # 1: true; 0: false
        ind = ind.astype('bool')
        ele_ind_list.append(ind)

    node_ind_list = []
    for i, shp in enumerate(shapes):
        print(f'shp {i+1} of {len(shapes)}')
        poly_xy = np.array(shp.points).T
        ind = inside_polygon(np.c_[grid.x, grid.y], poly_xy[0], poly_xy[1])  # 1: true; 0: false
        ind = ind.astype('bool')
        node_ind_list.append(ind)

    return [ele_ind_list, node_ind_list]
