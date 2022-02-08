from pylib import schism_grid, inside_polygon  # from ZG's pylib: pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple pylibs4schism==0.1.10
import matplotlib.pyplot as plt
import numpy as np
from time import time
from schism_py_pre_post.Geometry.inpoly import 
import os


if __name__ == "__main__":
    shapefile_name = f'/sciclone/schism10/feiye/Coastal_Act/Pumps/16a/levee_4_pump_polys.shp'
    maxelev = schism_grid('/sciclone/schism10/feiye/ICOGS/RUN21a/PostP/Maxelev/elev_max_stack121.1-121.24.gr3')
    hgrid = schism_grid('/sciclone/schism10/feiye/ICOGS/RUN21a/hgrid.gr3')

    # read pump capacities
    [ele_idx_list, gd] = find(
        # order matters, latter prevails
        shapefile_name=shapefile_name,
        hgrid_name=f"{w_dir}/hgrid.utm.gr3",
    )
    gd.compute_area()

    pass
