import numpy as np
from scipy import spatial
from schism_py_pre_post.Grid.SMS import find_pts_in_shpfiles
from pylib import schism_grid, sms2grd
import copy


if __name__ == '__main__':
    # -------------------------------------------------------------------------
    wdir = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/Hgrid/Shapefiles/Combined_load_bathy/'
    original_grid_dir = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/Hgrid/Shapefiles/Combined_load_bathy/'

    # thalweg
    thalweg_grd_name = f'{wdir}/merged_thalwegs_redist100m.2dm'
    # thalweg buffer
    thalweg_buffer_poly_shpfname = f'{wdir}/stream_redist_100m_buf_34m.shp'
    # full grid with bathymetry loaded
    full_grid_name = f'{original_grid_dir}/hgrid.utm.26918.gr3'
    # -------------------------------------------------------------------------

    full_gd = schism_grid(full_grid_name)
    thalweg_gd = sms2grd(thalweg_grd_name)

    fullgrid_thalweg_idx = find_pts_in_shpfiles([thalweg_buffer_poly_shpfname], pts=np.c_[full_gd.x, full_gd.y])

    thalweg_gd_pts = np.c_[thalweg_gd.x, thalweg_gd.y]
    fullgrid_thalweg_pts = np.c_[full_gd.x[fullgrid_thalweg_idx], full_gd.y[fullgrid_thalweg_idx]]

    mindist_thalweg_idx = spatial.cKDTree(thalweg_gd_pts).query(fullgrid_thalweg_pts)[1]
    full_gd.dp[fullgrid_thalweg_idx] = np.maximum(full_gd.dp[fullgrid_thalweg_idx], thalweg_gd.dp[mindist_thalweg_idx])

    print('writing test outputs ...\n')
    test_gd = copy.deepcopy(full_gd)
    test_gd.dp[:] = 0
    test_gd.dp[fullgrid_thalweg_idx] = 1
    test_gd.save(f'{wdir}/test.gr3')

    print('writing hgrid in utm ...\n')
    full_gd.save(f'{wdir}/hgrid.utm.26918.thalweg.gr3')
    full_gd.grd2sms(fname=f'{wdir}/hgrid.utm.26918.thalweg.2dm')

    print('writing hgrid in lon\lat ...\n')
    full_gd.x, full_gd.y = full_gd.proj(prj0='epsg:26918', prj1='epsg:4326')
    full_gd.save(f'{wdir}/hgrid.thalweg.ll')

    # print('writing hgrid in cpp ...\n')
    # full_gd.x, full_gd.y = full_gd.proj(prj0='epsg:4326', prj1='cpp', lon0=-77.07, lat0=24)
    # full_gd.save(f'{wdir}/hgrid.thalweg.cpp.gr3')
    # full_gd.write_hgrid(f'{wdir}/hgrid.thalweg.cpp')

    pass
