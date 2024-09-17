"""Test pybind11 module for SCHISM hgrid and vgrid file I/O."""


from pylib_experimental.schism_file import cread_schism_hgrid
from pylib import schism_grid
import numpy as np


gd_fname = ('/sciclone/schism10/feiye/STOFS3D-v6/'
            'Inputs/v6.1.static_inputs_snapshot20231105/hgrid.gr3')
gd_with_bnd = cread_schism_hgrid(gd_fname)

gd_fname = '/sciclone/schism10/feiye/STOFS3D-v6/Inputs/hgrid_no_boundary.gr3'
gd_without_bnd = cread_schism_hgrid(gd_fname)
gd_without_bnd.compute_bnd()

print('Done')
