import numpy as np
import os
from pylib import schism_grid, read_schism_bpfile, inside_polygon
import pandas as pd
from scipy.interpolate import griddata, Rbf
import netCDF4 as nc
import schism_py_pre_post
# import schism_py_pre_post.idw as idw
from schism_py_pre_post.Geometry.interp import inverse_distance_weighting


def gen_subregion_ic_stofs3d(wdir=None, obsdir=None, hycom_TS_file=None, date_str='2000-01-01'):
    # copy datafiles
    mydir = os.path.dirname(schism_py_pre_post.__file__)
    os.system(f'cp {mydir}/Datafiles/ecgc_sub_grid.reg {wdir}')
    os.system(f'cp {mydir}/Datafiles/ecgc_shoreline_sal.txt {wdir}')

    var_dict = {
        'tem': {'interp_method': 2, 'f_ecgc_shoreline': None},
        'sal': {'interp_method': 1, 'f_ecgc_shoreline': 'ecgc_shoreline_sal.txt'},
    }

    for var in ['tem', 'sal']:
        # ----------------------------------------------------------------------------------
        # -----------------inputs--------------------------
        # ----------------------------------------------------------------------------------
        # wdir = '/sciclone/schism10/feiye/From_Nabi/RUN02/Hotstart_v1/'
        interp_method = var_dict[var]['interp_method']  # recommended: 1 sal; 2 tem
        # USGS data, generated by download_usgs_with_api-----------------
        f_usgs_coastal = f'{obsdir}/mean_{var}_xyz_{date_str}'

        # Manually defined coastal shoreline, which will take HYCOM values,
        # The purpose is to improve the interpolation along the shoreline
        # Not need for temperature because the horizontal gradient is small compared to salinity
        if var_dict[var]['f_ecgc_shoreline'] is not None:
            f_ecgc_shoreline = wdir + var_dict[var]['f_ecgc_shoreline']
        else:
            f_ecgc_shoreline = None

        # HYCOM file, only used for interpolating values on f_ecgc_shoreline
        # hycom_TS_file = '/sciclone/schism10/feiye/Coastal_Act/Hot/13c/TS_1.nc'

        # Define sub region within which {var} is to be modified/x_ma
        f_grid_in = wdir + 'hgrid.ll'
        # this region is manually made, which covers the entire watershed area
        f_sub_region = wdir + 'ecgc_sub_grid.reg'

        # outputs:
        # same as f_grid_in, but dp will be interpolated from f_usgs_coastal and f_ecgc_shoreline within f_sub_region
        f_grid_out = wdir + 'ecgc_coastal_' + var + '.gr3'

        # ----------------------------------------------------------------------------------
        # -----------------read region--------------------------
        # ----------------------------------------------------------------------------------
        region = read_schism_bpfile(f_sub_region, fmt=1)

        # ----------------------------------------------------------------------------------
        # -----------------load txt--------------------------
        # ----------------------------------------------------------------------------------
        # usgs coastal watershed
        usgs_coastal = np.array(pd.read_csv(f_usgs_coastal, delimiter=' ', header=None).values[:, :-1]).astype(float)

        # shoreline
        if f_ecgc_shoreline is None:
            ecgc_shoreline = None
        else:
            ecgc_shoreline = np.loadtxt(f_ecgc_shoreline)

        # ----------------------------------------------------------------------------------
        # -----------------load netcdf--------------------------
        # ----------------------------------------------------------------------------------
        if f_ecgc_shoreline is not None:
            ds = nc.Dataset(hycom_TS_file)
            lon = ds['xlon']
            lat = ds['ylat']
            if var == 'sal':
                val = ds['salinity']
            elif var == 'tem':
                val = ds['temperature']

            [xx, yy] = np.meshgrid(lon, lat)
            m_mask = np.ma.getmask(val[0, 0, :, :])
            x_mask = np.ma.masked_array(xx, mask=m_mask)
            y_mask = np.ma.masked_array(yy, mask=m_mask)
            for i, _ in enumerate(ecgc_shoreline):
                dist = np.sqrt((x_mask - ecgc_shoreline[i, 0])**2 + (y_mask - ecgc_shoreline[i, 1])**2)
                ind = np.unravel_index(np.argmin(dist, axis=None), dist.shape)
                if ecgc_shoreline[i, 2] > 1e-6:
                    ecgc_shoreline[i, 2] = val[0, 0, ind[0], ind[1]]
            xyz = np.append(usgs_coastal, ecgc_shoreline, axis=0)
        else:
            xyz = usgs_coastal

        # ensure quality
        xyz_min = -99999
        xyz_max = 99999
        if var == 'sal':
            xyz_min = 0
            xyz_max = 42
        elif var == 'tem':
            xyz_min = -5
            xyz_max = 42
        box = [-150, -50, 0, 50]
        # remove invalid values
        idx = np.where((xyz[:, 2] >= xyz_min) & (xyz[:, 2] <= xyz_max))[0]
        xyz = xyz[idx, :]
        # remove invalid lon/lat
        idx = np.where((xyz[:, 0] >= box[0]) & (xyz[:, 0] <= box[1]))[0]
        xyz = xyz[idx, :]
        np.savetxt(f'{wdir}/ecgc_shoreline+coastal_{var}.txt', xyz, delimiter=' ')

        # ----------------------------------------------------------------------------------
        # -----------------load hgrid--------------------------
        # ----------------------------------------------------------------------------------
        gd = schism_grid(f_grid_in)
        idx = (inside_polygon(np.c_[gd.x, gd.y], region.x, region.y) == 1)
        gd_x = gd.x[idx]
        gd_y = gd.y[idx]

        # gd = Hgrid("/sciclone/schism10/feiye/ICOGS/RUN06f/Hot_and_3Dth_for_06h/ecgc.gr3")
        # gd.get_nodes()
        # gd_x = gd.node_np[:, 1]
        # gd_y = gd.node_np[:, 2]

        # ----------------------------------------------------------------------------------
        # -----------------interpolation--------------------------
        # ----------------------------------------------------------------------------------
        if interp_method == 0:  # Rbf interpolator
            interpolator = Rbf(xyz[:, 0], xyz[:, 1], xyz[:, 2], function='linear')
            z_interp = interpolator(gd_x, gd_y)
            z_interp = np.maximum(z_interp, min(xyz[:, 2]))
            z_interp = np.minimum(z_interp, max(xyz[:, 2]))
        elif interp_method == 1:  # linear
            z_interp = griddata(xyz[:, 0:2], xyz[:, 2], (gd_x, gd_y), method='linear')
        elif interp_method == 2:  # inverse-distance weighted
            # idw_tree = idw.tree(xyz[:, 0:2], xyz[:, 2])  # see github: paulbrodersen/inverse_distance_weighting
            # z_interp = idw_tree(np.transpose(np.vstack([gd_x, gd_y])))
            z_interp = inverse_distance_weighting(
                X=xyz[:, 0:2], val=xyz[:, 2],
                Xq=np.transpose(np.vstack([gd_x, gd_y])), n_nei=6
            )
        else:
            raise Exception(f'unknown interp_method {interp_method}')
        nan_idx = np.where(np.isnan(z_interp))

        # dealing with nan
        z_interp_nearest = griddata(xyz[:, 0:2], xyz[:, 2], (gd_x, gd_y), method='nearest')
        z_interp[nan_idx] = z_interp_nearest[nan_idx]

        gd.dp[:] = -9999
        gd.dp[idx] = z_interp
        gd.write_hgrid(f_grid_out)
