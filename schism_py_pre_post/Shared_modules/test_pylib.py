"""test pylib functions"""

from pylib import read_schism_output
from pylib_experimental.schism_file import TimeHistory
import numpy as np

data = read_schism_output(
    run='/sciclone/schism10/feiye/STOFS3D-v8/R03b/',
    varname=['elevation'],
    xyz='/sciclone/schism10/feiye/STOFS3D-v8/BPfiles/USGS_station_LA_repositioned.bp',
    stacks=np.arange(1, 36),
)

th = TimeHistory(data_array=np.c_[data.time, data.elevation.T])
th.writer('/sciclone/schism10/feiye/STOFS3D-v8/O03b/elevation.USGS_station_LA_repositioned.dat')

print("done!")
