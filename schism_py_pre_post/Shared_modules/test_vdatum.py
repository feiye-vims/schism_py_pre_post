from coastalmodeling_vdatum import vdatum, _path
from pylib import schism_grid

PATH = '/sciclone/schism10/Hgrid_projects/GEOTIFFs/'

_path.NAVD88_G2018=f"{PATH}/us_noaa_g2018u0.tif"
_path.XGEOID20B=f"{PATH}/xGEOID20B.tif"
_path.MLLW_ITRF2020_2020=f"{PATH}/us_noaa_nos_MLLW-ITRF2020_2020.0_nwldatum_4.7.0_20240621_.tif"
_path.LMSL_ITRF2020_2020=f"{PATH}/us_noaa_nos_LMSL-ITRF2020_2020.0_nwldatum_4.7.0_20240621_.tif"

gd = schism_grid('/sciclone/schism10/feiye/STOFS3D-v8/I09/hgrid.gr3')
x,y,z = vdatum.convert("navd88", "xgeoid20b", gd.y, gd.x, gd.z, online=False, epoch=None)
x, y, z = vdatum.convert("navd88", "xgeoid20b", 29.87809000, -91.31890000 , 0, online=False, epoch=None)


from vyperdatum.points import VyperPoints
import numpy as np


vp = VyperPoints(vdatum_directory='/sciclone/schism10/Hgrid_projects/DEMs/vdatum/vdatum_all_20220511')

vp = VyperPoints()
x = np.array([-76.19698, -76.194, -76.198])
y = np.array([37.1299, 37.1399, 37.1499])
z = np.array([10.5, 11.0, 11.5])

vp.transform_points((3631, 'mllw'), (6319, 'mllw'), x, y, z=z)
pass