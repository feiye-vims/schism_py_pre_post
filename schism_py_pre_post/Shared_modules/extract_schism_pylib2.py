"""extract schism output using pylib"""

import numpy as np
from pylib import read_schism_output

# --------------- inputs ----------------
# USGS_repositioned_v49 for v8, whole domain
# USGS_station_LA_repositioned_nontidal_v43_paper.bp for v8, LA, nontidal
# USGS_station_LA_repositioned_nontidal_v7p1_paper.bp, LA, nontidal

# RUN_DIR = '/sciclone/schism10/feiye/STOFS3D-v8/RUN16_v6/'
# output_files = [
#     '/sciclone/schism10/feiye/STOFS3D-v8/O16_v6/elevation.USGS_repositioned_v49.dat'
# ]
# bpfile = [
#     '/sciclone/schism10/feiye/STOFS3D-v8/BPfiles/USGS_repositioned_v49.bp'
# ]

# RUN_DIR = '/sciclone/schism10/feiye/STOFS3D-v8/R13_v7/'
# output_files = [
#     '/sciclone/schism10/feiye/STOFS3D-v8/O13_v7/elevation.USGS_station_LA_repositioned_nontidal.dat'
# ]
# bpfile = [
#     '/sciclone/schism10/feiye/STOFS3D-v8/BPfiles/USGS_station_LA_repositioned_nontidal_v43.bp'
# ]

# RUN_DIR = '/sciclone/schism10/feiye/STOFS3D-v8/R13p_v7/'
# output_files = [
#     '/sciclone/schism10/feiye/STOFS3D-v8/O13p_v7/outputs/elevation.coops_test.dat'
# ]
# bpfile = [
#     '/sciclone/schism10/feiye/STOFS3D-v8/BPfiles/coops_test.bp'
# ]

RUN_DIR = '/sciclone/schism10/feiye/STOFS3D-v8/R09k4/'
output_files = [
    '/sciclone/schism10/feiye/STOFS3D-v8/O09k4/elevation.USGS_station_LA_repositioned_nontidal_v43_paper.dat',
    # '/sciclone/schism10/feiye/STOFS3D-v8/O09f4/elevation.USGS_station_LA_tidal.dat',
]
bpfile = [
    '/sciclone/schism10/feiye/STOFS3D-v8/BPfiles/USGS_station_LA_repositioned_nontidal_v43_paper.bp',
    # '/sciclone/schism10/feiye/STOFS3D-v8/BPfiles/USGS_station_LA_tidal.bp',
]
# ---------------------------------------

start_stack = 1
end_stack = 55
for bpfile, output_file in zip(bpfile, output_files):
    data = read_schism_output(run=RUN_DIR, varname=['elevation'], xyz=bpfile, stacks=np.arange(start_stack, end_stack+1))

    np_savetxt_args = {'fmt': '%.4f', 'delimiter': ' ', 'newline': '\n'}
    np.savetxt(output_file, np.c_[data.time, data.elevation.T], **np_savetxt_args)

print("done!")
