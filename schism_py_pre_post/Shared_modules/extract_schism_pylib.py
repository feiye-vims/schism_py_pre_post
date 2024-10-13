"""extract schism output using pylib"""

import numpy as np
from pylib import read_schism_output
from pylib_experimental.schism_file import TimeHistory

RUN_DIR = '/sciclone/schism10/feiye/STOFS3D-v8/R09c/'
output_files = [
    '/sciclone/schism10/feiye/STOFS3D-v8/O09c/elevation.USGS_station_LA_repositioned_nontidal.dat',
    '/sciclone/schism10/feiye/STOFS3D-v8/O09c/elevation.USGS_station_LA_repositioned.dat'
]
bpfile = [
    '/sciclone/schism10/feiye/STOFS3D-v8/BPfiles/USGS_station_LA_repositioned_nontidal_v43.bp',
    '/sciclone/schism10/feiye/STOFS3D-v8/BPfiles/USGS_station_LA_repositioned_v43.bp'
]

for bpfile, output_file in zip(bpfile, output_files):
    data = read_schism_output(run=RUN_DIR, varname=['elevation'], xyz=bpfile, stacks=np.arange(1, 36))
    th = TimeHistory(data_array=np.c_[data.time, data.elevation.T])

    # from matplotlib import pyplot as plt
    # plt.plot(th.time, th.data)
    # plt.show()

    th.writer(output_file)

print("done!")
