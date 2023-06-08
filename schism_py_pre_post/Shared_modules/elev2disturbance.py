from schism_py_pre_post.Grid.Hgrid_extended import read_schism_hgrid_cached
from pylib import grd2sms
import copy

wdir = '/sciclone/schism10/feiye/Coastal_Act/RUN13h/'
gd = read_schism_hgrid_cached(f'{wdir}/hgrid.ll')
elevmax = read_schism_hgrid_cached(f'{wdir}/elev_max_stack53-60.gr3')

upland = gd.dp < 0
disturbance = copy.deepcopy(elevmax)
disturbance.dp[upland] += gd.dp[upland]

grd2sms(disturbance, f'{wdir}/disturbance_stack53-60.2dm')
pass
