import numpy as np
import os
import schism_py_pre_post
from schism_py_pre_post.Shared_modules.hotstart_proc import Hotstart
from schism_py_pre_post.Download.download_usgs_with_api import get_usgs_obs_for_stofs3d
from schism_py_pre_post.Download.download_cbp_with_api import get_cbp_obs_for_stofs3d
from schism_py_pre_post.Shared_modules.gen_subregion_ic2 import gen_subregion_ic_stofs3d
from schism_py_pre_post.Grid.Grid_geometry import find_points_in_polyshp
from pylib import schism_grid
from pathlib import Path


# input section
hotstart_date_str = '2023-06-01'
wdir = '/sciclone/schism10/feiye/STOFS3D-v6/Inputs/I_20230601_new/Hotstart/'
griddir = wdir
output_obs_dir = f'{wdir}/Obs/'
hycom_TS_file = f'{wdir}/TS_1.nc'
hycom_hot_file = f'{wdir}/hotstart.nc.hycom'
city_shapefile_names = ["city_polys_from_v10_lonlat.shp"]
my_hot_file = f'{wdir}/hotstart.nc'
# end input section

# copy datafiles
mydir = os.path.dirname(schism_py_pre_post.__file__)
for shp in city_shapefile_names:
    shp_basename = Path(shp).stem
    os.system(f'cp {mydir}/Datafiles/{shp_basename}.* {wdir}')

# download coastal obs from usgs
get_usgs_obs_for_stofs3d(outdir=output_obs_dir, start_date_str=hotstart_date_str)

# download coastal obs from CBP
get_cbp_obs_for_stofs3d(outdir=output_obs_dir, sample_time=hotstart_date_str)

# interpolate obs onto model grid
gen_subregion_ic_stofs3d(wdir=wdir, obsdir=output_obs_dir, hycom_TS_file=hycom_TS_file, date_str=hotstart_date_str)

# make a copy of the hycom-based hotstart.nc
if os.path.exists(my_hot_file):
    os.system(f"rm {my_hot_file}")
os.system(f"cp {hycom_hot_file} {my_hot_file}")

# tweak coastal values based on obs
my_hot = Hotstart(
    grid_info=griddir,
    hot_file=my_hot_file
)

for i, var in enumerate(['tem', 'sal']):
    hg = schism_grid(f'{wdir}/ecgc_coastal_{var}.gr3')  # this file is from get*_obs_for_stofs3d
    idx = hg.dp > -9998
    for k in range(my_hot.grid.vgrid.nvrt):
        my_hot.tr_nd.val[idx, k, i] = hg.dp[idx]

# set salinity to 0 on higher grounds
rat = np.maximum(np.minimum(1.0, (my_hot.grid.hgrid.dp + 3.0) / 3.0), 0.0)  # linearly varying from 0 to 3 m
my_hot.tr_nd.val[:, :, 1] *= np.transpose(np.tile(rat, (my_hot.grid.vgrid.nvrt, 1)))
my_hot.trnd_propogate()  # propogate trnd values to trnd0 and tr_el

# set initial elevation: 0 in the ocean, just below ground on higher grounds and in cities
h0 = 0.1
above_NAVD88_0 = my_hot.grid.hgrid.dp < 0
shapefile_names = [f'{wdir}/{x}' for x in city_shapefile_names]
in_city = find_points_in_polyshp(pt_xy=np.c_[hg.x, hg.y], shapefile_names=shapefile_names)

ind = np.logical_or(above_NAVD88_0, in_city)

elev_ic = np.zeros(my_hot.grid.hgrid.dp.shape)
elev_ic[ind] = - my_hot.grid.hgrid.dp[ind] - h0
my_hot.eta2.val[:] = elev_ic

# write
my_hot.writer(fname=my_hot_file)
