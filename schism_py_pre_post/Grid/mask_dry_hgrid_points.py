from schism_py_pre_post.Grid.Hgrid_ported import read_schism_hgrid_cached
from pylib import schism_grid, sms2grd, grd2sms


wdir = '/sciclone/schism10/feiye/STOFS3D-v6/Outputs/O13k/PostP/MaxElev/'
maxelev_name = 'elev_max_stack27-34'

gd = read_schism_hgrid_cached(f'{wdir}/hgrid.ll')
maxelev = read_schism_hgrid_cached(f'{wdir}/{maxelev_name}.gr3')

dry = maxelev.dp + gd.dp < 0.0
maxelev.dp[dry] = -9999.0
maxelev.save(f'{wdir}/{maxelev_name}.masked.gr3')
grd2sms(maxelev, f'{wdir}/{maxelev_name}.masked.2dm')
