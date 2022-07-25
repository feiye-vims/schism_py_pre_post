# %%
from unittest.main import main
from unittest.util import sorted_list_difference
from pylib import schism_grid, sms2grd, grd2sms
import os
import numpy as np

# %%  
def find_large_small_dp():
    gd = sms2grd('/sciclone/schism10/feiye/STOFS3D-v4/Shared_with_NOAA/Mask_lbnd/maxelev.mask.2dm')
    sorted_dp = np.sort(gd.dp)
    sorted_idx = np.argsort(gd.dp)
    valid = np.squeeze(np.argwhere(~np.isnan(sorted_dp)))
    n = 30
    print(np.c_[sorted_dp[valid[-n:]], sorted_idx[valid[-n:]], gd.x[sorted_idx[valid[-n:]]], gd.y[sorted_idx[valid[-n:]]]])
    with open('/sciclone/schism10/feiye/STOFS3D-v4/Shared_with_NOAA/Mask_lbnd/large_elev.txt', 'w') as file:  
        file.write('x y z id\n')
        np.savetxt(file, np.c_[gd.x[sorted_idx[valid[-n:]]], gd.y[sorted_idx[valid[-n:]]], sorted_dp[valid[-n:]], sorted_idx[valid[-n:]]])
    with open('/sciclone/schism10/feiye/STOFS3D-v4/Shared_with_NOAA/Mask_lbnd/small_elev.txt', 'w') as file:  
        file.write('x y z id\n')
        np.savetxt(file, np.c_[gd.x[sorted_idx[valid[:n]]], gd.y[sorted_idx[valid[:n]]], sorted_dp[valid[:n]], sorted_idx[valid[:n]]])


def proc_hgrid():
    file_2dm = f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/Hgrid/Shapefiles/Combined_1.2.2/v5_1.2.2.2dm'
    dirname = os.path.dirname(file_2dm)

    gd = sms2grd(file_2dm)
    utm_x, utm_y = [gd.x, gd.y]
    gd.compute_area()
    gd.compute_ctr()
    sorted_area = np.sort(gd.area)
    sorted_idx = np.argsort(gd.area)
    print(np.c_[sorted_area[:20], sorted_idx[:20], gd.xctr[sorted_idx[:20]], gd.yctr[sorted_idx[:20]]])

    gd.x, gd.y = gd.proj(prj0='epsg:26918', prj1='epsg:4326')
    gd.save(f'{dirname}/hgrid.ll')
    pass

    # load bathymetry
    # os.chdir(dirname)
    # os.system('mpirun ./pload_depth.py')

    os.system(f'mv {dirname}/hgrid.ll.new {dirname}/hgrid.ll')
    gd_ll_loaded = schism_grid(f'{dirname}/hgrid.ll')

    # gd.x, gd.y = gd_ll_loaded.proj(prj0='epsg:4326', prj1='cpp', lon0=-77.07, lat0=24)
    # gd.dp = gd_ll_loaded.dp
    # gd.write_hgrid(f'{dirname}/hgrid.cpp')

    gd.x = utm_x
    gd.y = utm_y
    gd.dp = gd_ll_loaded.dp
    gd.write_hgrid(f'{dirname}/hgrid.utm.26918')
    os.chdir(dirname)
    os.system('ln -s hgrid.utm.26918 hgrid.utm.26918.gr3')

    grd2sms(gd, file_2dm)
    pass

if __name__ == "__main__":
    proc_hgrid()
