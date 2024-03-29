from schism_py_pre_post.Geometry.inpoly import find_node_in_shpfiles
from pylib import schism_grid, proj
import numpy as np
import os
import copy


def set_additional_dp_v11_91(gd_ll=None, gd_dem=None, wdir='./'):
    # v11.91, set upstream mississippi river levee height to 25 m, to prevent overtopping near the source
    #         set Bonnet Carre Spill Way gate to 8.5 m
    #         revert levee height near Ostrica, LA to DEM values

    # Check levee info existence
    levee_info_dir = f'{wdir}/Levee_info/'
    if not os.path.exists(levee_info_dir):
        levee_tar_file = os.path.dirname(Datafiles.__file__)  + "/Levee_info.tar"
        my_tar = tarfile.open(levee_tar_file)
        my_tar.extractall(wdir)
        my_tar.close()


    if gd_ll is None:
        gd_ll = schism_grid(f'{wdir}/hgrid.ll')
    if gd_dem is None:
        gd_dem = schism_grid(f'{wdir}/hgrid.DEM_loaded.ll')
    
    gd_meters = copy.deepcopy(gd_ll)
    gd_meters.proj(prj0='epsg:4326', prj1='esri:102008')

    for set_depth, shapefile in zip(
        [-9, -20],
        [f'{levee_info_dir}/Additional_Polygons/BonnetCarre_102008.shp',
         f'{levee_info_dir}/Additional_Polygons/la_levee_center_line_upstream_missi_13m_buffer_102008.shp']
    ):
        i_inpoly = find_node_in_shpfiles(shapefile_names=[shapefile], gd=gd_meters)
        gd_ll.dp[i_inpoly] = np.minimum(gd_ll.dp[i_inpoly], set_depth)

    i_inpoly = find_node_in_shpfiles(shapefile_names=[f'{levee_info_dir}/Additional_Polygons/la_levee_center_line_Ostrica_buffer_102008.shp'], gd=gd_meters)
    gd_ll.dp[i_inpoly] = gd_dem.dp[i_inpoly]


    # os.system(f"cp {wdir}/hgrid.additional_dp.ll {wdir}/hgrid.ll")
    # proj(
    #     f'{wdir}/hgrid.ll', 0, 'epsg:4326',
    #     f'{wdir}/hgrid.utm.gr3', 0, 'epsg:26918',
    # )
    # os.system(f"cp {wdir}/hgrid.utm.gr3 {wdir}/hgrid.utm.26918")

    return gd_ll

if __name__ == "__main__":
    # needs hgrid.ll in wdir
    wdir = '/sciclone/schism10/feiye/STOFS3D-v6/Inputs/V6_mesh_from_JZ/Test_levee_heights/'
    gd_ll = set_additional_dp_v11_91(wdir=wdir)
    gd_ll.save(f'{wdir}/hgrid.additional_dp.ll')