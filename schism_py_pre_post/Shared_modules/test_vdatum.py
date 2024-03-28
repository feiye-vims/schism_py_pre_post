from vyperdatum.points import VyperPoints
import numpy as np


vp = VyperPoints(vdatum_directory='/sciclone/schism10/Hgrid_projects/DEMs/vdatum/vdatum_all_20220511')

vp = VyperPoints()
x = np.array([-76.19698, -76.194, -76.198])
y = np.array([37.1299, 37.1399, 37.1499])
z = np.array([10.5, 11.0, 11.5])

vp.transform_points((3631, 'mllw'), (6319, 'mllw'), x, y, z=z)
pass