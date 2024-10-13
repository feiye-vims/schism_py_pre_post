from pylib import inside_polygon
import shapefile
import numpy as np


def find_points_in_polyshp(pt_xy, shapefile_names):
    ind = np.zeros(pt_xy[:, 0].shape)
    for shapefile_name in shapefile_names:
        # shapefile # records = sf.records() # sf.shapeType # len(sf) # s = sf.shape(idx)
        sf = shapefile.Reader(shapefile_name)
        shapes = sf.shapes()

        for i, shp in enumerate(shapes):
            poly_xy = np.array(shp.points).T
            print(f'shape {i+1} of {len(shapes)}, {poly_xy[:, 0]}')
            ind += inside_polygon(pt_xy, poly_xy[0], poly_xy[1])  # 1: true; 0: false

    ind = ind.astype('bool')
    return ind


def sample():
    from pylib_experimental.schism_file import cread_schism_hgrid
    gd = cread_schism_hgrid('/sciclone/schism10/feiye/STOFS3D-v8/I09b/Drag/drag.gr3')
    idx = find_points_in_polyshp(
        np.c_[gd.x, gd.y], ['/sciclone/schism10/feiye/STOFS3D-v8/I09b/Drag/atcha_mouth.shp']
    )
    gd.dp[idx] *= 10  # scale z by 10
    gd.save('/sciclone/schism10/feiye/STOFS3D-v8/I09b/Drag/drag_tweaked.gr3')


if __name__ == '__main__':
    sample()