"""extract schism output using pylib"""

import numpy as np
from pylib import read_schism_output
from pylib_experimental.schism_file import TimeHistory

# --------------- inputs ----------------
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

RUN_DIR = '/sciclone/schism10/feiye/STOFS3D-v8/R13p_v7/'
output_files = [
    '/sciclone/schism10/feiye/STOFS3D-v8/O13p_v7/outputs/elevation.coops_test.dat'
]
bpfile = [
    '/sciclone/schism10/feiye/STOFS3D-v8/BPfiles/coops_test.bp'
]
# ---------------------------------------

start_stack = 1
end_stack = 396
for bpfile, output_file in zip(bpfile, output_files):
    data = read_schism_output(run=RUN_DIR, varname=['elevation'], xyz=bpfile, stacks=np.arange(start_stack, end_stack+1))
    th = TimeHistory(data_array=np.c_[data.time, data.elevation.T])

    # from matplotlib import pyplot as plt
    # plt.plot(th.time, th.data)
    # plt.show()

    np_savetxt_args = {'fmt': '%.4f', 'delimiter': ' ', 'newline': '\n'}
    np.savetxt(output_file, np.c_[th.time + start_stack-1 + 1/24, th.data], **np_savetxt_args)
    # th.writer(output_file)

print("done!")
