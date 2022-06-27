#!/usr/bin/env python3

import shapefile
# import matplotlib.pyplot as plt
from pylib import schism_grid, inside_polygon
import numpy as np
from time import time
import schism_py_pre_post.Shared_modules.inpoly_operations as inpoly_operations


def set_val_gr3_shpfile(shapefile_names: list = [""],
                        ops: list = [inpoly_operations.set_val], op_kwargs=None,
                        op_background=None, op_bg_kwarg=None,
                        hgrid_name: str = "", hgrid_xy_name: str = "", outfilename: str = ""):
    '''
    Change the z value in gr3 within all polygons of a shapefile.
    Also takes a list of shapefiles and different values can be assigned to each shapefile

    i_operation:
        3 types of operations can be specified:
            -1: gd.dp=min(gd.dp, val)
            0: gd.dp=val
            1: gd.dp=max(gd.dp, val)
        If necessary, cjreate additional user-defined operations in "def operation"

    hgrid_name:
        input hgrid to be modified

    hgrid_xy_name:
        Optionally, a different hgrid can be specified, whose x,y will be used for
        the "in_polygon" operation but not in the output grid
    '''
    # hgrid
    t = time()
    gd = schism_grid(hgrid_name)
    if (hgrid_xy_name is None) or (hgrid_name == hgrid_xy_name):
        gd_xy = gd
    else:
        gd_xy = schism_grid(hgrid_xy_name)
    print(f'Reading hgrid took {time()-t} seconds')

    # set default value outside polygons
    dp = gd.dp.copy()
    if op_background is not None:
        if op_bg_kwarg is None:
            dp[:] = op_background(gd.dp)
        else:
            dp[:] = op_background(gd.dp, None, **op_bg_kwarg)

    # modify hgrid within regions specified in shapefiles
    for op_func, op_kwargs, shapefile_name in zip(ops, op_kwargs, shapefile_names):
        t = time()
        # shapefile # records = sf.records() # sf.shapeType # len(sf) # s = sf.shape(idx)
        sf = shapefile.Reader(shapefile_name)
        shapes = sf.shapes()

        ind = np.zeros(dp.shape)
        for shp in shapes:
            poly_xy = np.array(shp.points).T
            ind += inside_polygon(np.c_[gd_xy.x, gd_xy.y], poly_xy[0], poly_xy[1])  # 1: true; 0: false
        ind = ind.astype('bool')
        dp[ind] = op_func(gd.dp[ind], ind, **op_kwargs)

        print(f'Processing {shapefile_name} took {time()-t} seconds')

    gd.dp = dp.copy()
    # optional: plot
    # gd.plot_grid(method=0, fmt=1)
    # plt.show()

    gd.write_hgrid(fname=outfilename)
    pass


if __name__ == "__main__":
    set_val_gr3_shpfile(
        # order matters, latter prevails
        shapefile_names=["/sciclone/schism10/feiye/ICOGS/Ida04b/Shapiro_mod/shp/carri_islands_ll_buf_25km.shp",
                         "/sciclone/schism10/feiye/ICOGS/Ida04b/Shapiro_mod/shp/carri_islands_ll_buf_10km.shp",
                         "/sciclone/schism10/feiye/ICOGS/Ida04b/Shapiro_mod/shp/carri_islands_ll.shp"],
        ops=[inpoly_operations.add_val,
             inpoly_operations.add_val,
             inpoly_operations.add_val],
        op_kwargs=[{'val_spec': 0.05, 'val_thres': 0.5},
                   {'val_spec': 0.1, 'val_thres': 0.5},
                   {'val_spec': 0.3, 'val_thres': 0.5}],
        op_background=None, op_bg_kwarg=None,
        hgrid_name="/sciclone/schism10/feiye/ICOGS/Ida04b/shapiro.gr3",
        hgrid_xy_name="/sciclone/schism10/feiye/ICOGS/Ida04b/hgrid.pkl",  # used for determining "in_polygon", should have the same projection as shpfiles
        outfilename='/sciclone/schism10/feiye/ICOGS/Ida04b/Shapiro_mod/shapiro.mod.gr3'
    )
