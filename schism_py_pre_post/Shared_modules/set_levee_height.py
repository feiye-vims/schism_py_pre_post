import numpy as np
from pylib import schism_grid
from schism_py_pre_post import Datafiles
from schism_py_pre_post.Geometry.inpoly import find_node_in_shpfiles
import copy
import os
import shutil
import tarfile


def set_constant_levee_height(gd=None, wdir='./'):
    # set constant levee heights for levees not in NLD.

    # Check levee info existence
    levee_info_dir = f'{wdir}/Levee_info/'
    if not os.path.exists(levee_info_dir):
        levee_tar_file = os.path.dirname(Datafiles.__file__)  + "/Levee_info.tar"
        my_tar = tarfile.open(levee_tar_file)
        my_tar.extractall(wdir)
        my_tar.close()

    if gd is None:
        gd = schism_grid(f'{wdir}/hgrid.ll')
    
    gd_meters = copy.deepcopy(gd)
    gd_meters.proj(prj0='epsg:4326', prj1='esri:102008')

    for set_depth, shapefile in zip(
        [-9, -9], [f'{levee_info_dir}/Polygons/additional_levee_poly_102008.shp', f"{levee_info_dir}/Polygons/HH_Dam_buffer_13m_102008.shp"]
    ):
        i_inpoly = find_node_in_shpfiles(shapefile_names=[shapefile], gd=gd_meters)
        gd.dp[i_inpoly] = np.minimum(gd.dp[i_inpoly], set_depth)

    return gd

if __name__ == "__main__":
    # needs hgrid.ll in wdir
    # needs Levee_info/, copy it into the wdir
    wdir = '/sciclone/schism10/feiye/STOFS3D-v6/Inputs/V6_mesh_from_JZ/Test_levee_heights/'
    gd = set_constant_levee_height(wdir=wdir)
    gd.save(f'{wdir}/hgrid.constant_levee_height.ll')
