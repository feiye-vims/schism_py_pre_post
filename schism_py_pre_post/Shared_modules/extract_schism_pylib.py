"""extract schism output using pylib"""

import numpy as np
from pylib import read_schism_output
from pylib_experimental.schism_file import TimeHistory

run_dir = '/sciclone/schism10/feiye/STOFS3D-v8/R07b/'
output_file = '/sciclone/schism10/feiye/STOFS3D-v8/O07b/elevation.USGS_station_LA_repositioned_nontidal'
bpfile = '/sciclone/schism10/feiye/STOFS3D-v8/BPfiles/USGS_station_LA_repositioned_nontidal.bp'

data = read_schism_output(run=run_dir, varname=['elevation'], xyz=bpfile, stacks=np.arange(1, 36))

th = TimeHistory(data_array=np.c_[data.time, data.elevation.T])

# from matplotlib import pyplot as plt
# plt.plot(th.time, th.data)
# plt.show()

th.writer(output_file)

print("done!")
