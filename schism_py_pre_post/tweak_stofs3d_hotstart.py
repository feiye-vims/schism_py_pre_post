import numpy as np
import os
from schism_py_pre_post.hotstart_proc import Hotstart
from schism_py_pre_post.Download.download_usgs_with_api import get_usgs_obs_for_stofs3d
from schism_py_pre_post.Download.download_cbp_with_api import get_cbp_obs_for_stofs3d
from schism_py_pre_post.gen_subregion_ic2 import gen_subregion_ic_stofs3d
from pylib import schism_grid


# input section
wdir = '/sciclone/schism10/feiye/From_Nabi/RUN02/Hotstart_v1/'
griddir = '/sciclone/schism10/feiye/From_Nabi/RUN02/Hotstart_v1/'
output_obs_dir = '/sciclone/schism10/feiye/From_Nabi/RUN02/Hotstart_v1/Obs/'
hycom_hot_file = '/sciclone/schism10/feiye/From_Nabi/RUN02/Hotstart_v1/hotstart.nc.0'
my_hot_file = '/sciclone/schism10/feiye/From_Nabi/RUN02/Hotstart_v1/hotstart.nc'
hotstart_date_str = '2015-09-18'
# end input section

mypath = os.path.dirname(os.path.realpath(__file__))
os.system(f'cp {mypath}/ecgc_shoreline_sal.txt {wdir}')
os.system(f'cp {mypath}/ecgc_sub_grid.reg {wdir}')

# download coastal obs from usgs
# get_usgs_obs_for_stofs3d(outdir=output_obs_dir)

# download coastal obs from CBP
# get_cbp_obs_for_stofs3d(outdir=output_obs_dir)

# make a copy of the hycom-based hotstart.nc
if os.path.exists(my_hot_file):
    os.system(f"rm {my_hot_file}")
os.system(f"cp {hycom_hot_file} {my_hot_file}")

# tweak coastal values
my_hot = Hotstart(
    grid_info=griddir,
    hot_file=my_hot_file
)

for i, var in enumerate(['tem', 'sal']):
    hg = schism_grid(f'{wdir}/ecgc_coastal_{var}.gr3')  # this file is from get*_obs_for_stofs3d
    idx = hg.dp > -9998
    for k in range(my_hot.grid.vgrid.nvrt):
        my_hot.tr_nd.val[idx, k, i] = hg.dp[idx]

# set salinity to 0 for higher grounds
rat=np.maximum(np.minimum(1.0, (my_hot.grid.hgrid.dp+3.0)/3.0), 0.0)
my_hot.tr_nd.val[:, :, 1] *= np.transpose(np.tile(rat, (my_hot.grid.vgrid.nvrt, 1)))

# write
my_hot.writer(fname=my_hot_file)
