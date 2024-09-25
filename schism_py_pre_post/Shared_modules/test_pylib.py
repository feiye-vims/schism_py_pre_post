"""test pylib functions"""

# %%
from pylib import read
bp1 = read('/sciclone/schism10/feiye/STOFS3D-v8/BPfiles/USGS_station_LA_repositioned_nontidal.bp')
bp2 = read('/sciclone/schism10/feiye/STOFS3D-v8/BPfiles/USGS_station_LA_repositioned_v43.bp')

for i, station in enumerate(bp1.station):
    idx = bp2.station == station
    if sum(idx) == 0:
        raise ValueError(f"Station {station} not found in bp2.")
    elif sum(idx) > 1:
        raise ValueError(f"Station {station} found more than once in bp2.")
    bp1.x[i] = bp2.x[idx]
    bp1.y[i] = bp2.y[idx]
    bp1.z[i] = bp2.z[idx]

bp1.save('/sciclone/schism10/feiye/STOFS3D-v8/BPfiles/USGS_station_LA_repositioned_nontidal_v43.bp')

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
