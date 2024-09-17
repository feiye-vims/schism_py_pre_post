"""test pylib functions"""

# %%
from pylib import convert_schism_source
convert_schism_source(run='/sciclone/schism10/feiye/STOFS3D-v7/Shadow_fcst/',)

# %%
from pylib import read_schism_output
from pylib_experimental.schism_file import TimeHistory
import numpy as np

data = read_schism_output(
    run='/sciclone/schism10/feiye/STOFS3D-v7/Runs/R15a/',
    varname=['elevation'],
    xyz='/sciclone/schism10/feiye/STOFS3D-v8/BPfiles/USGS_station_LA_repositioned_nontidal.bp',
    stacks=np.arange(1, 36),
)

th = TimeHistory(data_array=np.c_[data.time, data.elevation.T])

# from matplotlib import pyplot as plt
# plt.plot(th.time, th.data)
# plt.show()

th.writer('/sciclone/schism10/feiye/STOFS3D-v7/Outputs/O15a/elevation.USGS_station_LA_repositioned_nontidal.dat')

print("done!")

# %%
