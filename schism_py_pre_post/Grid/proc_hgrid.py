# %%
from unittest.util import sorted_list_difference
from pylib import schism_grid, sms2grd, grd2sms
import os
import numpy as np

# %%
file_2dm = f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/Hgrid/Shapefiles/LA1.0/LA1.0.2dm'
dirname = os.path.dirname(file_2dm)

# %%
gd = sms2grd(file_2dm)
utm_x, utm_y = [gd.x, gd.y]
gd.compute_area()
gd.compute_ctr()
sorted_area = np.sort(gd.area)
sorted_idx = np.argsort(gd.area)
print(np.c_[sorted_area[:20], sorted_idx[:20], gd.xctr[sorted_idx[:20]], gd.yctr[sorted_idx[:20]]])


# %%
gd.x, gd.y = gd.proj(prj0='epsg:26918', prj1='epsg:4326')
gd.save(f'{dirname}/hgrid.ll')
pass

# %% load bathymetry
# os.chdir(dirname)
# os.system('mpirun ./pload_depth.py')

# %%
os.system(f'mv {dirname}/hgrid.ll.new {dirname}/hgrid.ll')
gd_ll_loaded = schism_grid(f'{dirname}/hgrid.ll')

# gd.x, gd.y = gd_ll_loaded.proj(prj0='epsg:4326', prj1='cpp', lon0=-77.07, lat0=24)
# gd.dp = gd_ll_loaded.dp
# gd.write_hgrid(f'{dirname}/hgrid.cpp')

gd.x = utm_x
gd.y = utm_y
gd.dp = gd_ll_loaded.dp
gd.write_hgrid(f'{dirname}/hgrid.utm.26918')

grd2sms(gd, file_2dm)
pass
