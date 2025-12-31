"""
Check source/sink files by plotting cm/h at each element
"""

import numpy as np
from pylib_experimental.schism_file import SourceSink, cread_schism_hgrid as read_schism_hgrid

wdir = '/sciclone/schism10/feiye/STOFS3D-v8/I15i/Source_sink/Manual_combine_source_sink/'

ss = SourceSink.from_files(wdir)

hg = read_schism_hgrid(f'{wdir}/hgrid.gr3')
hg.proj(prj0='EPSG:4326', prj1='esri:102008')
hg.compute_area(); hg.compute_ctr()

sink_rate = ss.vsink.data / hg.area[ss.sink_eles - 1] * 3600 * 100  # m^3/s / m^2 * 3600 * 100 = cm/h

np.savetxt(
    f'{wdir}/sink_rate_cm_per_h.txt',
    np.c_[hg.xctr[ss.sink_eles - 1], hg.yctr[ss.sink_eles - 1], sink_rate.mean(axis=0)], fmt='%f %f %f',
    header='x y sink_rate_cm_per_h'
)

print('Done')
