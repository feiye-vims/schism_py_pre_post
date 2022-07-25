import xarray
import numpy as np
from pylib import schism_grid, grd2sms
from schism_py_pre_post.Geometry.inpoly import find_node_in_shpfiles
import os

mask_value = -9999

hg = schism_grid('/sciclone/schism10/feiye/STOFS3D-v4/Inputs/v4_20220715_update/hgrid.ll')
hg.compute_bnd()
hg.compute_node_ball()
# hg.save('/sciclone/schism10/feiye/STOFS3D-v4/Inputs/v4_20220715_update/hgrid.pkl')

lbnd_nd_idx = np.hstack(hg.ilbn)
ilbnd_nd = np.zeros((hg.np), dtype=bool)
ilbnd_nd[lbnd_nd_idx] = True

lbnd_side_idx = np.sum(ilbnd_nd[hg.isidenode], axis=1).astype(bool)

ilbnd_ele = np.unique(np.reshape(hg.isdel[lbnd_side_idx], (-1, )))
ilbnd_ele = ilbnd_ele[ilbnd_ele>=0]

lbnd_nd_2_tier_idx = np.unique(np.reshape(hg.elnode[ilbnd_ele], (-1, )))
lbnd_nd_2_tier_idx = lbnd_nd_2_tier_idx[lbnd_nd_2_tier_idx>=0].astype(int)
ilbnd_nd_2_tier = np.zeros((hg.np), dtype=bool)
ilbnd_nd_2_tier[lbnd_nd_2_tier_idx] = True

# find high ground
i_high_ground = hg.dp < -19.0

# get additional mask points
i_additional = find_node_in_shpfiles(shapefile_names=['/sciclone/schism10/feiye/STOFS3D-v4/Shared_with_NOAA/Mask_lbnd/additional_masks.shp'], gd=hg)

imask = i_high_ground + ilbnd_nd_2_tier + i_additional

ds = xarray.open_dataset('/sciclone/schism10/feiye/STOFS3D-v4/Shared_with_NOAA/Mask_lbnd/lbnd_node_mask.nc')
ds['idmask'][:] = imask.astype(int)
ds.to_netcdf('lbnd_node_mask.nc')
os.system('mv lbnd_node_mask.nc /sciclone/schism10/feiye/STOFS3D-v4/Shared_with_NOAA/Mask_lbnd/lbnd_node_mask.nc')
ds.close()

# generate out2d_1.mask.nc using nco commands
# /sciclone/schism10/feiye/STOFS3D-v4/Shared_with_NOAA/Mask_lbnd/mask_lbn_nd.sh

maxelev = schism_grid('/sciclone/schism10/feiye/STOFS3D-v4/Shared_with_NOAA/Mask_lbnd/maxelev.gr3')
maxelev.dp[imask] = np.nan
grd2sms(grd=maxelev, sms='/sciclone/schism10/feiye/STOFS3D-v4/Shared_with_NOAA/Mask_lbnd/maxelev.mask.2dm')

ds = xarray.open_dataset('/sciclone/schism10/feiye/STOFS3D-v4/Shared_with_NOAA/Mask_lbnd/out2d_1.mask.nc')
hg.dp[:] = ds['elevation'][-1, :]
ds.close()
grd2sms(grd=hg, sms='/sciclone/schism10/feiye/STOFS3D-v4/Shared_with_NOAA/Mask_lbnd/out2d_1.mask.2dm')